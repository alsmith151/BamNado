use ahash::AHashSet as HashSet;
use anyhow::{Ok, Result};
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use log::{debug, error, info, log, warn};
use noodles::bam::record;
use noodles::bed::record::OptionalFields;
use noodles::core::Position;
use noodles::{bam, bed, sam};
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use sam::alignment::record::data::field::value::Value;
use std::collections::HashMap;
use std::io;
use std::io::Write;
use std::path::PathBuf;
use std::prelude::rust_2021::*;
use std::str;
use tempfile;

use crate::utils::progress_bar;
use crate::utils::CB;

type Iv = Interval<usize, u32>;

struct ChromosomeRange {
    chrom: String,
    start: u64,
    end: u64,
}

#[derive(Debug, Clone)]
struct BedGraphEntry {
    chrom: String,
    start: u64,
    end: u64,
    score: f64,
}

#[derive(Debug)]
struct ChromStats {
    chrom: String,
    length: u64,
    mapped: u64,
    unmapped: u64,
}

#[derive(Debug)]
struct BamStats {
    mapped: u64,
    unmapped: u64,
    stats: HashMap<String, ChromStats>,
}

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum NormalizationMethod {
    Raw,
    RPKM,
    CPM,
}

impl NormalizationMethod {
    pub fn from_str(s: &str) -> Result<NormalizationMethod> {
        match s.to_lowercase().as_str() {
            "raw" => Ok(NormalizationMethod::Raw),
            "rpkm" => Ok(NormalizationMethod::RPKM),
            "cpm" => Ok(NormalizationMethod::CPM),
            _ => {
                warn!("Unknown normalization method: {}. Defaulting to Raw", s);
                Ok(NormalizationMethod::Raw)
            }
        }
    }

    pub fn scale_factor(&self, scale_factor: f32, bin_size: u64, n_reads: u64) -> f64 {
        match self {
            Self::Raw => 1.0,
            Self::CPM => {
                let total_reads = n_reads as f64 / 1e6;
                let scale_factor = scale_factor as f64 * total_reads;
                scale_factor
            }
            Self::RPKM => {
                let total_reads = n_reads as f64 / 1e6;
                let tile_len_kb = (bin_size as f64 / 1000 as f64);
                let reads_per_tile = total_reads as f64 * tile_len_kb as f64;
                let reads_per_tile = 1 as f64 / reads_per_tile;
                let scale_factor = scale_factor as f64 * reads_per_tile;
                scale_factor
            }
        }
    }
}

fn get_bam_header(file_path: PathBuf) -> Result<sam::Header> {
    // Check that the file exists
    if !file_path.exists() {
        return Err(anyhow::Error::from(io::Error::new(
            io::ErrorKind::NotFound,
            format!("File not found: {}", file_path.display()),
        )));
    };

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(file_path.clone())
        .expect("Failed to open file");

    let header = match reader.read_header() {
        std::result::Result::Ok(header) => header,
        Err(e) => {
            debug!(
                "Failed to read header using noodels falling back to samtools: {}",
                e
            );

            let header_samtools = std::process::Command::new("samtools")
                .arg("view")
                .arg("-H")
                .arg(file_path.clone())
                .output()
                .expect("Failed to run samtools")
                .stdout;

            let header_str =
                String::from_utf8(header_samtools).expect("Failed to convert header to string");

            // Slight hack here for CellRanger BAM files that are missing the version info
            let header_string =
                header_str.replace("@HD\tSO:coordinate\n", "@HD\tVN:1.6\tSO:coordinate\n");
            let header_str = header_string.as_bytes();
            let mut reader = sam::io::Reader::new(header_str);
            let header = reader
                .read_header()
                .expect("Failed to read header with samtools");
            header
        }
    };
    Ok(header)
}

struct BamDetails {
    file_path: PathBuf,
    is_single_cell: bool,
}

impl BamDetails {
    pub fn new(file_path: PathBuf, is_single_cell: bool) -> Self {
        Self {
            file_path: file_path,
            is_single_cell: is_single_cell,
        }
    }

    pub fn stats(&self) -> Result<BamStats> {
        let reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(self.file_path.clone())
            .expect("Failed to open file");
        let index = reader.index();
        let header = get_bam_header(self.file_path.clone())?;

        let mut chrom_stats = HashMap::new();

        for ((reference_sequence_name_buf, reference_sequence), index_reference_sequence) in header
            .reference_sequences()
            .iter()
            .zip(index.reference_sequences())
        {
            let reference_sequence_name = str::from_utf8(reference_sequence_name_buf)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            let (mapped_record_count, unmapped_record_count) = index_reference_sequence
                .metadata()
                .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                .unwrap_or_default();

            let stats = ChromStats {
                chrom: reference_sequence_name.to_string(),
                length: usize::from(reference_sequence.length()) as u64,
                mapped: mapped_record_count,
                unmapped: unmapped_record_count,
            };

            chrom_stats.insert(reference_sequence_name.to_string(), stats);
        }

        let mut unmapped_record_count = index.unplaced_unmapped_record_count().unwrap_or_default();

        let mut mapped_record_count = 0;

        for (_, stats) in chrom_stats.iter() {
            mapped_record_count += stats.mapped;
            unmapped_record_count += stats.unmapped;
        }

        Ok(BamStats {
            mapped: mapped_record_count,
            unmapped: unmapped_record_count,
            stats: chrom_stats,
        })
    }

    pub fn scale_factor(
        &self,
        scale_factor: f32,
        bin_size: u64,
        method: NormalizationMethod,
    ) -> Result<f64> {
        // Calculate the scale factor used to normalize the coverage values
        let stats = self.stats()?;
        let scale_factor = method.scale_factor(scale_factor, bin_size, stats.mapped);
        Ok(scale_factor)
    }

    fn genome_chunk_length(&self, bin_size: u64) -> Result<u64> {
        // genomeLength = sum(bamHandles[0].lengths)
        // max_reads_per_bp = max([float(x) / genomeLength for x in mappedList])
        // genomeChunkLength = int(min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandles)))))
        // genomeChunkLength -= genomeChunkLength % tile_size

        let stats = self.stats()?;
        let genome_length = stats.stats.values().map(|x| x.length).sum::<u64>();
        let max_reads_per_bp = stats.mapped as f64 / genome_length as f64;
        let genome_chunk_length = f64::from(5e7).min(2e6 / (max_reads_per_bp));

        let correction = genome_chunk_length % bin_size as f64;
        let genome_chunk_length = genome_chunk_length - correction;

        Ok(genome_chunk_length as u64)
    }

    pub fn chromsizes(&self) -> Result<HashMap<String, u64>> {
        let stats = match self.stats() {
            std::result::Result::Ok(stats) => stats,
            Err(e) => {
                error!("Failed to get stats: {}", e);
                return Err(e);
            }
        };

        let chrom_sizes: HashMap<String, u64> = stats
            .stats
            .iter()
            .map(|(chrom, stats)| (chrom.clone(), stats.length))
            .collect();
        Ok(chrom_sizes)
    }

    pub fn write_chromsizes(&self, outfile: PathBuf) -> Result<()> {
        let chrom_sizes = self.chromsizes()?;

        info!("Writing chromosome sizes to file...");
        let mut writer = io::BufWriter::new(std::fs::File::create(outfile)?);

        for (chrom, length) in chrom_sizes.iter() {
            writeln!(writer, "{}\t{}", chrom, length)?;
        }
        Ok(())
    }
}

struct IntervalMaker<'a> {
    use_fragment: bool,
    is_single_cell: bool,
    details: &'a Option<SingleCellDetails>,
}

impl<'a> IntervalMaker<'a> {
    fn create_interval(&self, record: bam::Record, use_fragment: bool) -> Option<Iv> {
        let record = record;
        let is_paired = !record.flags().is_mate_unmapped();
        let is_read_1 = record.flags().is_first_segment();
        let is_reverse = record.flags().is_reverse_complemented();
        let is_mate_reverse = record.flags().is_mate_reverse_complemented();
        let is_proper_pair = record.flags().is_properly_segmented();

        // Two ways of generating the start and end of the interval:
        // 1. Using just the read start and end
        // 2. Using the fragment:
        //
        //     R1                      R2
        //    |----------------->     <-----------------|
        //    r1_start               r2_end
        //
        //    R2                      R1
        //    |------------------->    <-----------------|
        //    r2_start                  r1_end
        //
        let coords = match use_fragment {
            true => {
                // Check if the read is the first in the segment. Dont want to double count
                if is_read_1 {
                    // Check if the read is properly paired
                    if is_proper_pair {
                        // Check if the read is in the forward orientation
                        if !is_reverse {
                            let r1_start = record
                                .alignment_start()
                                .expect("Failed to get start")
                                .expect("Error with getting start")
                                .get();
                            let r1_length = record.sequence().len();
                            let r1_end = r1_start + r1_length;

                            let start = r1_end;
                            let end = start + record.template_length() as usize;
                            Some((start, end))
                        } else {
                            let r2_start = record
                                .mate_alignment_start()
                                .expect("Failed to get mate start")
                                .expect("Error with getting mate start")
                                .get();

                            let start = r2_start;
                            let end = start + record.template_length().abs() as usize;
                            Some((start, end))
                        }
                    } else {
                        // If the read is not properly paired, return None
                        None
                    }
                } else {
                    // If the read is not the first in the segment, return None
                    // This is to avoid double counting
                    None
                }
            }
            false => {
                // If not using the fragment, just use the read start and end
                let start = record
                    .alignment_start()
                    .expect("Failed to get start")
                    .expect("Error with getting start")
                    .get();

                // The end is the start + the length of the read sequence
                let end = record.sequence().len() + start;
                Some((start, end))
            }
        };

        // Create the interval
        let interval = match coords {
            Some((start, end)) => Some(Iv {
                start: start,
                stop: end,
                val: 0,
            }),
            None => None,
        };

        interval
    }

    fn create_interval_single_cell(
        &self,
        record: bam::Record,
        use_fragment: bool,
        details: &SingleCellDetails,
    ) -> Option<Iv> {
        // Get the barcode from the record
        let rec = record.clone();
        let name = rec.name().expect("Failed to get name");

        let barcode = match rec.data().get(&CB) {
            Some(b) => {
                let barcode = match b {
                    std::result::Result::Ok(Value::String(b)) => str::from_utf8(b)
                        .expect("Failed to get barcode")
                        .to_string(),
                    _ => {
                        warn!("Barcode is not a string in record {:?}", name);
                        return None;
                    }
                };
                barcode
            }
            None => {
                debug!("No barcode found in record {:?}", name);
                return None;
            }
        };

        // Check if the barcode is in the list of barcodes
        match details.barcodes.contains(&barcode) {
            true => self.create_interval(record, use_fragment),
            false => None,
        }
    }

    pub fn interval(&self, record: bam::Record) -> Option<Iv> {
        match (self.is_single_cell, self.details) {
            (true, Some(details)) => {
                self.create_interval_single_cell(record, self.use_fragment, details)
            }
            _ => self.create_interval(record, self.use_fragment),
        }
    }
}

#[derive(Debug, Clone)]
pub struct SingleCellDetails {
    barcodes: HashSet<String>,
    barcode_flag: [u8; 2],
}

impl SingleCellDetails {
    pub fn new(barcodes: HashSet<String>, barcode_flag: Option<[u8; 2]>) -> Self {
        let barcode_flag = match barcode_flag {
            Some(flag) => flag,
            None => CB,
        };
        Self {
            barcodes: barcodes,
            barcode_flag: barcode_flag,
        }
    }
}

pub struct BamPileup {
    file_path: PathBuf,
    bin_size: u64,
    norm_method: NormalizationMethod,
    use_fragment: bool,
    scale_factor: f32,
    single_cell_details: Option<SingleCellDetails>,
}

impl BamPileup {
    pub fn new(
        file_path: PathBuf,
        bin_size: u64,
        use_fragment: bool,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        single_cell_details: Option<SingleCellDetails>,
    ) -> Self {
        Self {
            file_path: file_path,
            bin_size: bin_size,
            norm_method: norm_method,
            use_fragment: use_fragment,
            scale_factor: scale_factor,
            single_cell_details: single_cell_details,
        }
    }

    fn get_genomic_chunks(&self, details: BamDetails) -> Result<Vec<ChromosomeRange>> {
        let genome_chunk_length = details.genome_chunk_length(self.bin_size)?;
        let chrom_stats = details.stats()?;
        let chrom_chunks: Vec<ChromosomeRange> = chrom_stats
            .stats
            .iter()
            .map(|(chrom, stats)| {
                let mut chunks = Vec::new();
                let mut start = 1;
                let end = stats.length;
                while start <= end {
                    // Corrected to include the last position in the range
                    let chunk_end = start + genome_chunk_length - 1; // Adjust to ensure the chunk covers exactly genome_chunk_length positions
                    let chunk_end = chunk_end.min(end); // Ensure we do not exceed the chromosome length
                    chunks.push(ChromosomeRange {
                        chrom: chrom.to_string(),
                        start: start,
                        end: chunk_end,
                    });
                    start = chunk_end + 1; // Corrected to start the next chunk right after the current chunk ends
                }
                chunks
            })
            .flatten()
            .collect();
        Ok(chrom_chunks)
    }

    fn pileup(&self, norm_method: NormalizationMethod) -> Result<Vec<BedGraphEntry>> {
        let bam_details =
            BamDetails::new(self.file_path.clone(), self.single_cell_details.is_some());
        let scale_factor =
            bam_details.scale_factor(self.scale_factor, self.bin_size, norm_method.clone())?;
        let chrom_chunks = self.get_genomic_chunks(bam_details)?;

        // Options for counting
        let chunks_len = chrom_chunks.len();
        let normalise = match norm_method {
            NormalizationMethod::Raw => false,
            _ => true,
        };
        let use_fragment = self.use_fragment;
        let is_single_cell = self.single_cell_details.is_some();
        let single_cell_details = self.single_cell_details.clone();
        let bedgraph: Vec<BedGraphEntry> = chrom_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                chunks_len as u64,
                "Performing pileup".to_string(),
            ))
            .map(|chrom_range| {
                // Open the file

                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(self.file_path.clone())
                    .expect("Failed to open file");
                let header = get_bam_header(self.file_path.clone()).expect("Failed to get header");

                // Extract the region
                let coords = format!(
                    "{}:{}-{}",
                    chrom_range.chrom, chrom_range.start, chrom_range.end
                );
                let region = coords.parse().expect("Failed to parse region");

                // Generate the bins for the region (bin_size)
                let starts = (chrom_range.start..=chrom_range.end).step_by(self.bin_size as usize); // Corrected to include the end in the range
                let ends = starts
                    .clone()
                    .skip(1)
                    .chain(std::iter::once(chrom_range.end + 1)); // Adjusted to ensure the last bin includes the end position
                let n_bins =
                    ((chrom_range.end - chrom_range.start) + self.bin_size - 1) / self.bin_size; // Adjusted to correctly calculate the number of bins when the range is not a perfect multiple of bin_size

                // Build a Lapper object from the reads in the region
                let records: Box<bam::io::reader::Query<noodles::bgzf::Reader<std::fs::File>>> =
                    reader
                        .query(&header, &region)
                        .map(Box::new)
                        .expect("Failed to query region");

                let intervals: Vec<Iv> = records
                    .map(|record| {
                        let interval = IntervalMaker {
                            use_fragment: use_fragment,
                            is_single_cell: is_single_cell,
                            details: &single_cell_details,
                        }
                        .interval(record.expect("Failed to get record"));
                        interval
                    })
                    .filter_map(|x| x)
                    .collect();

                let bam_lapper = Lapper::new(intervals);

                // Count overlaps with the bins
                // If the score is the same as the previous entry, extend the previous entry
                let mut bedgraph: Vec<BedGraphEntry> = Vec::with_capacity(n_bins as usize);
                for (start, end) in starts.zip(ends) {
                    let count = bam_lapper.count(start as usize, end as usize);
                    let score = match normalise {
                        true => count as f64 * scale_factor,
                        false => count as f64,
                    };

                    let previous_entry = bedgraph.last_mut();
                    match previous_entry {
                        Some(prev) => {
                            if prev.score == score {
                                prev.end = end as u64 - 1;
                                continue;
                            }
                        }
                        None => {}
                    }

                    bedgraph.push(BedGraphEntry {
                        chrom: chrom_range.chrom.clone(),
                        start: start as u64,
                        end: end as u64 - 1,
                        score: score,
                    });
                }

                bedgraph
            })
            .flatten()
            .collect();

        Ok(bedgraph)
    }

    fn _sorted_pileup(&self) -> Result<Vec<BedGraphEntry>> {
        let bedgraph = self
            .pileup(self.norm_method.clone())?
            .into_iter()
            .sorted_by_key(|x| (x.chrom.clone(), x.start));

        Ok(bedgraph.collect())
    }

    fn _normalised_pileup(&self) -> Result<Vec<BedGraphEntry>> {
        let bedgraph = self._sorted_pileup()?;

        let n_total_reads = bedgraph.par_iter().map(|x| x.score).sum::<f64>();
        let scale_factor =
            self.norm_method
                .scale_factor(self.scale_factor, self.bin_size, n_total_reads as u64);

        let bedgraph = bedgraph
            .into_par_iter()
            .map(|entry| {
                let score = entry.score * scale_factor;
                BedGraphEntry {
                    chrom: entry.chrom,
                    start: entry.start,
                    end: entry.end,
                    score: score,
                }
            })
            .collect();

        Ok(bedgraph)
    }

    pub fn write_bedgraph(&self, outfile: PathBuf) -> Result<()> {
        // Extract the pileup and sort by chrom and start
        log::info!("Piling up the reads and sorting the bedgraph...");
        let bedgraph = self._normalised_pileup()?;
        let n_records = bedgraph.len();

        // Write the bedgraph to a file
        log::info!("Writing the bedgraph to a file...");
        let mut writer =
            bed::io::Writer::new(std::io::BufWriter::new(std::fs::File::create(outfile)?));

        for entry in bedgraph.into_iter().progress_with(progress_bar(
            n_records as u64,
            "Writing BedGraph".to_string(),
        )) {
            let optional_fields = OptionalFields::from(vec![entry.score.to_string()]);
            let record = bed::Record::<3>::builder()
                .set_reference_sequence_name(entry.chrom)
                .set_start_position(
                    Position::new(entry.start as usize)
                        .expect(format!("Failed to set start position {}", entry.start).as_str()),
                )
                .set_end_position(
                    Position::new(entry.end as usize)
                        .expect(format!("Failed to set end position {}", entry.end).as_str()),
                )
                .set_optional_fields(optional_fields.clone())
                .build()?;
            writer.write_record(&record)?;
        }

        Ok(())
    }

    pub fn write_bigwig(&self, outfile: PathBuf, n_threads: usize) -> Result<()> {
        log::info!("Obtaining BAM details...for conversion to BigWig");

        // BAM details
        let bam_details =
            BamDetails::new(self.file_path.clone(), self.single_cell_details.is_some());

        info!("Got BAM details");

        // Chromosome sizes
        let tempfile_chromsizes_path = tempfile::NamedTempFile::new()?
            .into_temp_path()
            .to_path_buf();

        info!("Setting up the chromosome sizes...");
        bam_details.write_chromsizes(tempfile_chromsizes_path.clone())?;

        // Write the bedgraph to a temporary file
        let tempfile_bdg_path = tempfile::NamedTempFile::new()?
            .into_temp_path()
            .to_path_buf();

        info!("Setting up the bedgraph...");
        self.write_bedgraph(tempfile_bdg_path.clone())?;

        // Use the bedGraphToBigWig tool to convert the bedgraph to bigwig as long as it is in the PATH
        // Convert the bedgraph to bigwig

        log::info!("Converting bedgraph to bigwig...");

        let bw_command_output = std::process::Command::new("bedGraphToBigWig")
            .arg(tempfile_bdg_path)
            .arg(tempfile_chromsizes_path)
            .arg(outfile.clone())
            .output()
            .expect("Failed to convert bedgraph to bigwig");

        // Check if the command was successful
        if !bw_command_output.status.success() {
            log::error!("Failed to convert bedgraph to bigwig");
            log::error!("{}", String::from_utf8_lossy(&bw_command_output.stderr));
            std::process::exit(1);
        } else {
            let message = format!("BigWig file written to {}", outfile.display());
            info!("{}", message);
        }

        Ok(())
    }
}
