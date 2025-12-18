//! # BAM Utilities Module
//!
//! This module contains a collection of utility functions and data structures for working with
//! BAM files. It includes tools for:
//! *   Parsing and handling BAM headers and indexes.
//! *   Managing cell barcodes and strand information.
//! *   Calculating basic statistics from BAM files.
//! *   Helper types for genomic intervals and coordinate systems.

use ahash::{HashMap, HashSet};
use anyhow::{Context, Result};
use bio_types::annot::contig::Contig;
use bio_types::strand::ReqStrand;
use log::debug;

use bio_types::strand::Strand as BioStrand;
use noodles::core::{Position, Region};
use noodles::{bam, bed, sam};
use polars::prelude::*;
use rust_lapper::Interval;
use rust_lapper::Lapper;
use std::ops::Bound;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// The cell barcode tag (CB).
pub const CB: [u8; 2] = [b'C', b'B'];
/// Type alias for an interval with a `u32` value.
pub type Iv = Interval<usize, u32>;

/// Represents the strand of a genomic feature.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum Strand {
    /// Forward strand (+).
    Forward,
    /// Reverse strand (-).
    Reverse,
    /// Both strands (unstranded).
    Both,
}
impl FromStr for Strand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "forward" => Ok(Strand::Forward),
            "reverse" => Ok(Strand::Reverse),
            "both" => Ok(Strand::Both),
            _ => Err(format!("Invalid strand value: {s}")),
        }
    }
}
impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "forward"),
            Strand::Reverse => write!(f, "reverse"),
            Strand::Both => write!(f, "both"),
        }
    }
}

impl From<Strand> for BioStrand {
    fn from(val: Strand) -> Self {
        match val {
            Strand::Forward => BioStrand::Forward,
            Strand::Reverse => BioStrand::Reverse,
            Strand::Both => BioStrand::Unknown,
        }
    }
}

/// Creates a progress bar with a standard style.
///
/// # Arguments
///
/// * `length` - The total number of items to process.
/// * `message` - The message to display next to the progress bar.
pub fn progress_bar(length: u64, message: String) -> indicatif::ProgressBar {
    // Progress bar
    let progress_style = indicatif::ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}",
    )
    .unwrap()
    .progress_chars("##-");

    indicatif::ProgressBar::new(length)
        .with_style(progress_style)
        .with_message(message)
}

/// A collection of unique cell barcodes.
pub struct CellBarcodes {
    barcodes: HashSet<String>,
}

impl Default for CellBarcodes {
    fn default() -> Self {
        Self::new()
    }
}

impl CellBarcodes {
    /// Creates a new empty `CellBarcodes` collection.
    pub fn new() -> Self {
        Self {
            barcodes: HashSet::default(),
        }
    }

    /// Returns a copy of the set of barcodes.
    pub fn barcodes(&self) -> HashSet<String> {
        self.barcodes.clone()
    }

    /// Returns `true` if the collection is empty.
    pub fn is_empty(&self) -> bool {
        self.barcodes.is_empty()
    }

    /// Loads cell barcodes from a CSV file.
    ///
    /// The CSV file must have a header with a column named "barcode".
    pub fn from_csv(file_path: &PathBuf) -> Result<Self> {
        let path = Path::new(file_path).to_path_buf();

        let df = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path))?
            .finish()?;
        let mut barcodes = HashSet::default();

        for barcode in df.column("barcode").unwrap().str().unwrap() {
            let barcode = barcode.unwrap().to_string();
            barcodes.insert(barcode);
        }

        println!("Number of barcodes: {}", barcodes.len());
        println!(
            "First 10 barcodes: {:?}",
            barcodes.iter().take(10).collect::<Vec<_>>()
        );

        Ok(Self { barcodes })
    }
}

/// A collection of sets of cell barcodes, where each set represents a group.
pub struct CellBarcodesMulti {
    barcodes: Vec<HashSet<String>>,
}

impl Default for CellBarcodesMulti {
    fn default() -> Self {
        Self::new()
    }
}

impl CellBarcodesMulti {
    /// Creates a new empty `CellBarcodesMulti` collection.
    pub fn new() -> Self {
        Self {
            barcodes: Vec::new(),
        }
    }

    /// Returns a copy of the vector of barcode sets.
    pub fn barcodes(&self) -> Vec<HashSet<String>> {
        self.barcodes.clone()
    }

    /// Returns `true` if the collection is empty.
    pub fn is_empty(&self) -> bool {
        self.barcodes.is_empty()
    }

    /// Loads multiple sets of barcodes from a CSV file.
    ///
    /// Each column in the CSV file is treated as a separate set of barcodes.
    pub fn from_csv(file_path: &PathBuf) -> Result<Self> {
        let path = Path::new(file_path).to_path_buf();

        let df = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path))?
            .finish()?;

        let mut barcodes = Vec::new();

        for column in df.iter() {
            let barcode = column.str().unwrap();

            let bc = barcode
                .iter()
                .filter_map(|x| x.map(|x| x.to_string()))
                .collect::<HashSet<String>>();

            barcodes.push(bc);
        }

        Ok(Self { barcodes })
    }
}

#[derive(Debug)]
#[allow(dead_code)]
struct ChromosomeStats {
    chrom: String,
    length: u64,
    mapped: u64,
    unmapped: u64,
}

/// Reads the header of a BAM file.
///
/// If the `noodles` crate fails to read the header, it falls back to `samtools`.
/// This is useful for handling BAM files with slightly non-standard headers (e.g., from CellRanger).
pub fn bam_header(file_path: PathBuf) -> Result<sam::Header> {
    // Read the header of a BAM file
    // If the noodles crate fails to read the header, we fall back to samtools
    // This is a bit of a hack but it works for now
    // The noodles crate is more strict about the header format than samtools

    // Check that the file exists
    if !file_path.exists() {
        return Err(anyhow::Error::from(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("File not found: {}", file_path.display()),
        )));
    };

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(file_path.clone())
        .expect("Failed to open file");

    let header = match reader.read_header() {
        std::result::Result::Ok(header) => header,
        Err(e) => {
            debug!("Failed to read header using noodels falling back to samtools: {e}");

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
            reader
                .read_header()
                .expect("Failed to read header with samtools")
        }
    };
    Ok(header)
}

/// Stores statistics and metadata about a BAM file.
pub struct BamStats {
    // BAM file stats

    // File path
    #[allow(dead_code)]
    file_path: PathBuf,

    // Header
    header: sam::Header,

    // Contigs present in the BAM file
    #[allow(dead_code)]
    contigs: Vec<Contig<String, ReqStrand>>,

    // Chromosome stats
    chrom_stats: HashMap<String, ChromosomeStats>,

    // Total number of reads
    n_reads: u64,
    // Number of reads mapped to the genome
    n_mapped: u64,
    // Number of reads unmapped to the genome
    n_unmapped: u64,
}

impl BamStats {
    /// Creates a new `BamStats` instance from a BAM file path.
    ///
    /// This reads the header and index to calculate statistics.
    pub fn new(file_path: PathBuf) -> Result<Self> {
        // Read the header of the BAM file
        let header = bam_header(file_path.clone())?;

        // Get the contigs from the header
        let contigs = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| {
                let name = name.to_string();
                let length = map.length().get();
                Contig::new(name, 1, length, ReqStrand::Forward)
            })
            .collect();

        // Get the index of the BAM file
        let bam_reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(file_path.clone())
            .expect("Failed to open file");
        let index = bam_reader.index();

        // Get the chromosome stats
        let mut chrom_stats = HashMap::default();

        for ((reference_sequence_name_buf, reference_sequence), index_reference_sequence) in header
            .reference_sequences()
            .iter()
            .zip(index.reference_sequences())
        {
            let reference_sequence_name = std::str::from_utf8(reference_sequence_name_buf)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

            let (mapped_record_count, unmapped_record_count) = index_reference_sequence
                .metadata()
                .map(|m| (m.mapped_record_count(), m.unmapped_record_count()))
                .unwrap_or_default();

            let stats = ChromosomeStats {
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

        let n_reads = mapped_record_count + unmapped_record_count;

        Ok(Self {
            file_path,
            header,
            contigs,
            chrom_stats,
            n_reads,
            n_mapped: mapped_record_count,
            n_unmapped: unmapped_record_count,
        })
    }

    /// Estimates the optimal genome chunk length for parallel processing.
    ///
    /// The chunk length is calculated based on the genome size, number of mapped reads,
    /// and the desired bin size.
    pub fn estimate_genome_chunk_length(&self, bin_size: u64) -> Result<u64> {
        // genomeLength = sum(bamHandles[0].lengths)
        // max_reads_per_bp = max([float(x) / genomeLength for x in mappedList])
        // genomeChunkLength = int(min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandles)))))
        // genomeChunkLength -= genomeChunkLength % tile_size

        let stats = &self.chrom_stats;
        let genome_length = stats.values().map(|x| x.length).sum::<u64>();
        let max_reads_per_bp = self.n_mapped as f64 / genome_length as f64;
        let genome_chunk_length = 5e6_f64.min(2e6 / (max_reads_per_bp));

        let correction = genome_chunk_length % bin_size as f64;
        let genome_chunk_length = genome_chunk_length - correction;

        Ok(genome_chunk_length as u64)
    }

    /// Splits the genome into chunks for parallel processing.
    pub fn genome_chunks(&self, bin_size: u64) -> Result<Vec<Region>> {
        let genome_chunk_length = self.estimate_genome_chunk_length(bin_size)?;
        let chrom_chunks = self
            .chrom_stats
            .iter()
            .flat_map(|(chrom, stats)| {
                let mut chunks = Vec::new();
                let mut chunk_start = 1;
                let chrom_end = stats.length;

                while chunk_start <= chrom_end {
                    // Corrected to include the last position in the range
                    let chunk_end = chunk_start + genome_chunk_length - 1; // Adjust to ensure the chunk covers exactly genome_chunk_length positions
                    let chunk_end = chunk_end.min(chrom_end); // Ensure we do not exceed the chromosome length

                    let start = Position::try_from(chunk_start as usize).unwrap();
                    let end = Position::try_from(chunk_end as usize).unwrap();

                    let region = Region::new(&*chrom.clone(), start..=end);
                    chunks.push(region);
                    chunk_start = chunk_end + 1; // Corrected to start the next chunk right after the current chunk ends
                }
                chunks
            })
            .collect();

        Ok(chrom_chunks)
    }

    /// Splits a specific chromosome into chunks.
    pub fn chromosome_chunks(&self, chrom: &str, bin_size: u64) -> Result<Vec<Region>> {
        let genome_chunk_length = self.estimate_genome_chunk_length(bin_size)?;
        let stats = self
            .chrom_stats
            .get(chrom)
            .ok_or_else(|| anyhow::Error::msg(format!("Chromosome {chrom} not found")))?;

        let mut chunks = Vec::new();
        let mut chunk_start = 1;
        let chrom_end = stats.length;

        while chunk_start <= chrom_end {
            // Corrected to include the last position in the range
            let chunk_end = chunk_start + genome_chunk_length - 1; // Adjust to ensure the chunk covers exactly genome_chunk_length positions
            let chunk_end = chunk_end.min(chrom_end); // Ensure we do not exceed the chromosome length

            let start = Position::try_from(chunk_start as usize).unwrap();
            let end = Position::try_from(chunk_end as usize).unwrap();

            let region = Region::new(chrom, start..=end);
            chunks.push(region);
            chunk_start = chunk_end + 1; // Corrected to start the next chunk right after the current chunk ends
        }

        Ok(chunks)
    }

    /// Returns a mapping from chromosome ID to chromosome name.
    pub fn chromosome_id_to_chromosome_name_mapping(&self) -> HashMap<usize, String> {
        let mut ref_id_mapping = HashMap::default();
        for (i, (name, _)) in self.header.reference_sequences().iter().enumerate() {
            let name = std::str::from_utf8(name)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
                .unwrap()
                .to_string();
            ref_id_mapping.insert(i, name);
        }
        ref_id_mapping
    }

    /// Returns a mapping from chromosome name to chromosome ID.
    pub fn chromosome_name_to_id_mapping(&self) -> Result<HashMap<String, usize>> {
        let mut mapping = HashMap::default();
        let id_to_name = self.chromosome_id_to_chromosome_name_mapping();

        // Invert the mapping
        for (id, name) in id_to_name.iter() {
            mapping.insert(name.clone(), *id);
        }

        Ok(mapping)
    }

    /// Returns a mapping from chromosome ID to chromosome length.
    pub fn chromsizes_ref_id(&self) -> Result<HashMap<usize, u64>> {
        let mut chromsizes = HashMap::default();
        for (i, (_, map)) in self.header.reference_sequences().iter().enumerate() {
            let length = map.length().get();
            chromsizes.insert(i, length as u64);
        }
        Ok(chromsizes)
    }

    /// Returns a mapping from chromosome name to chromosome length.
    pub fn chromsizes_ref_name(&self) -> Result<HashMap<String, u64>> {
        let mut chromsizes = HashMap::default();
        for (name, map) in self.header.reference_sequences() {
            let length = map.length().get();
            let name = std::str::from_utf8(name)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
                .unwrap()
                .to_string();
            chromsizes.insert(name, length as u64);
        }
        Ok(chromsizes)
    }

    /// Writes chromosome sizes to a TSV file.
    pub fn write_chromsizes(&self, outfile: PathBuf) -> Result<()> {
        let chrom_sizes = self.chromsizes_ref_name()?;
        // Convert the chromosome sizes to a DataFrame
        let mut df = DataFrame::new(vec![
            Column::new(
                "chrom".into(),
                chrom_sizes.keys().cloned().collect::<Vec<String>>(),
            ),
            Column::new(
                "length".into(),
                chrom_sizes.values().cloned().collect::<Vec<u64>>(),
            ),
        ])?;

        // Sort the DataFrame by chromosome name
        df.sort_in_place(["chrom"], SortMultipleOptions::default())?;

        let mut file = std::fs::File::create(outfile)?;
        CsvWriter::new(&mut file)
            .include_header(false)
            .with_separator(b'\t')
            .finish(&mut df)?;
        Ok(())
    }

    /// Returns the number of mapped reads.
    pub fn n_mapped(&self) -> u64 {
        self.n_mapped
    }

    /// Returns the number of unmapped reads.
    pub fn n_unmapped(&self) -> u64 {
        self.n_unmapped
    }

    /// Returns the total number of reads.
    pub fn n_total_reads(&self) -> u64 {
        self.n_reads
    }
}

/// Supported output file types.
pub enum FileType {
    /// BedGraph format.
    Bedgraph,
    /// BigWig format.
    Bigwig,
    /// Tab-Separated Values.
    TSV,
}

impl std::str::FromStr for FileType {
    type Err = String;

    fn from_str(s: &str) -> Result<FileType, String> {
        match s.to_lowercase().as_str() {
            "bedgraph" => Ok(FileType::Bedgraph),
            "bdg" => Ok(FileType::Bedgraph),
            "bigwig" => Ok(FileType::Bigwig),
            "bw" => Ok(FileType::Bigwig),
            "tsv" => Ok(FileType::TSV),
            "txt" => Ok(FileType::TSV),
            _ => Err(format!("Unknown file type: {s}")),
        }
    }
}

/// Converts a list of regions to a `Lapper` interval tree for efficient overlap queries.
pub fn regions_to_lapper(regions: Vec<Region>) -> Result<HashMap<String, Lapper<usize, u32>>> {
    let mut lapper: HashMap<String, Lapper<usize, u32>> = HashMap::default();
    let mut intervals: HashMap<String, Vec<Iv>> = HashMap::default();

    for reg in regions {
        let chrom = reg.name().to_string();
        let start = reg.start().map(|x| x.get());
        let end = reg.end().map(|x| x.get());

        let start = match start {
            Bound::Included(start) => start,
            Bound::Excluded(start) => start + 1,
            _ => 0,
        };

        let end = match end {
            Bound::Included(end) => end,
            Bound::Excluded(end) => end - 1,
            _ => 0,
        };

        let iv = Iv {
            start,
            stop: end,
            val: 0,
        };
        intervals.entry(chrom.clone()).or_default().push(iv);
    }

    for (chrom, ivs) in intervals.iter() {
        let lap = Lapper::new(ivs.to_vec());
        lapper.insert(chrom.clone(), lap);
    }

    Ok(lapper)
}

/// Reads a BED file and converts it to a `Lapper` interval tree.
pub fn bed_to_lapper(bed: PathBuf) -> Result<HashMap<String, Lapper<usize, u32>>> {
    // Read the bed file and convert it to a lapper
    let reader = std::fs::File::open(bed)?;
    let buf_reader = std::io::BufReader::new(reader);
    let mut bed_reader = bed::io::Reader::<4, _>::new(buf_reader);
    let mut record = bed::Record::default();
    let mut intervals: HashMap<String, Vec<Iv>> = HashMap::default();
    let mut lapper: HashMap<String, Lapper<usize, u32>> = HashMap::default();

    while bed_reader.read_record(&mut record)? != 0 {
        let chrom = record.reference_sequence_name().to_string();
        let start = record.feature_start()?;
        let end = record
            .feature_end()
            .context("Failed to get feature end")??;

        let iv = Iv {
            start: start.get(),
            stop: end.get(),
            val: 0,
        };
        intervals.entry(chrom.clone()).or_default().push(iv);
    }

    for (chrom, ivs) in intervals.iter() {
        let lap = Lapper::new(ivs.to_vec());
        lapper.insert(chrom.clone(), lap);
    }
    // Check if the lapper is empty
    if lapper.is_empty() {
        return Err(anyhow::Error::msg("Lapper is empty"));
    }

    Ok(lapper)
}

/// Converts the chromosome names in a `Lapper` to chromosome IDs based on a BAM file's header.
///
/// This is useful for when we want to use the lapper with the BAM file.
pub fn convert_lapper_chrom_names_to_ids(
    lapper: HashMap<String, Lapper<usize, u32>>,
    bam_stats: &BamStats,
) -> Result<HashMap<usize, Lapper<usize, u32>>> {
    // Convert the chromosome names in the lapper to the chromosome IDs in the BAM file
    let chrom_id_mapping = bam_stats.chromosome_name_to_id_mapping()?;
    let mut lapper_chrom_id: HashMap<usize, Lapper<usize, u32>> = HashMap::default();
    for (chrom, lap) in lapper.iter() {
        let chrom_id = chrom_id_mapping.get(chrom).ok_or_else(|| {
            anyhow::Error::msg(format!("Chromosome {chrom} not found in BAM file"))
        })?;
        lapper_chrom_id.insert(*chrom_id, lap.clone());
    }
    Ok(lapper_chrom_id)
}

/// Reads the header of a BAM file.
///
/// If the `noodles` crate fails to read the header, it falls back to `samtools`.
/// This is useful for handling BAM files with slightly non-standard headers (e.g., from CellRanger).
pub fn get_bam_header<P>(file_path: P) -> Result<sam::Header>
where
    P: AsRef<Path>,
{
    let file_path = file_path.as_ref();

    // Check that the file exists
    if !file_path.exists() {
        return Err(anyhow::Error::from(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("File not found: {}", file_path.display()),
        )));
    };

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(file_path)
        .expect("Failed to open file");

    let header = match reader.read_header() {
        std::result::Result::Ok(header) => header,
        Err(e) => {
            debug!("Failed to read header using noodels falling back to samtools: {e}");

            let header_samtools = std::process::Command::new("samtools")
                .arg("view")
                .arg("-H")
                .arg(file_path)
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
            reader
                .read_header()
                .expect("Failed to read header with samtools")
        }
    };
    Ok(header)
}
