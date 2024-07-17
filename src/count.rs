use ahash::HashMap;
use ahash::HashMapExt;
use anyhow::{Result};
use bio_types::annot::loc::Loc;
use indicatif::ParallelProgressIterator;
use noodles::core::region::Region;
use noodles::core::Position;
use noodles::{bam, bed, core, sam};
use polars::prelude::*;
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::path::PathBuf;
use tempfile;
use log::{info, error, warn};
use std::fmt::Display;

use crate::filter::BamReadFilter;
use crate::intervals::IntervalMaker;
use crate::utils::{progress_bar, BamStats, Iv, NormalizationMethod};

struct BedgraphRecord {
    start: u64,
    end: u64,
    score: f32,
}

pub struct BamPileup {
    // The input BAM file
    file_path: PathBuf,
    
    // Genomic bin size
    bin_size: u64,

    // Normalization method to use for the pileup
    norm_method: NormalizationMethod,

    // Scaling factor for the pileup
    scale_factor: f32,

    // Use the fragment or the read for counting
    use_fragment: bool,

    // Filter used for the pileup
    filter: BamReadFilter,
}

impl Display for BamPileup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Pileup Settings: \n")?;
        write!(f, "BAM file: {}\n", self.file_path.display())?;
        write!(f, "Bin size: {}\n", self.bin_size)?;
        write!(f, "Normalizing using method: {:?}\n", self.norm_method)?;
        write!(f, "Scaling factor: {}\n", self.scale_factor)?;
        write!(f, "Using fragment for counting: {}\n", self.use_fragment)?;
        write!(f, "Filtering using: \n{}\n", self.filter)?;
        Ok(())

    }
}

impl BamPileup {
    pub fn new(
        file_path: PathBuf,
        bin_size: u64,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        use_fragment: bool,
        filter: BamReadFilter,
    ) -> Self {
        Self {
            file_path,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
            filter,
        }
    }

    // HashMap<String, Vec<bed::Record::<3>>>
    fn pileup(&self) -> Result<HashMap<String, Vec<Iv>>> {
        // Determine the chunks over which to iterate
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let genomic_chunks = bam_stats.genome_chunks(self.bin_size)?;

        // Set up required variables
        let chromsizes_refid = bam_stats
            .chromsizes_ref_id()
            .expect("Error getting chromsizes");

        // Iterate over the genomic chunks and pileup the reads
        let n_total_chunks = genomic_chunks.len();

        // Log important information
        // info!("Pileup settings:");
        // info!("BAM file: {}", self.file_path.display());
        // info!("Bin size: {}", self.bin_size);
        // info!("Performing pileup over {} genomic chunks", n_total_chunks);
        // info!("Normalizing using method: {:?}", self.norm_method);
        // info!("Scaling factor: {}", self.scale_factor);
        // info!("Using fragment for counting: {}", self.use_fragment);
        // info!("Filtering using: {}", self.filter);
        info!("{}", self.to_string());
        

        let pileup = genomic_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                n_total_chunks as u64,
                "Performing pileup".to_string(),
            ))
            .map(|region| {
                // Open the BAM file
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(self.file_path.clone())
                    .expect("Error opening BAM file");
                // Extract the header
                let header = reader.read_header().expect("Error reading BAM header");

                // Fetch the reads in the region
                let records = reader
                    .query(&header, &region)
                    .expect("Error querying BAM file");

                // Make intervals from the reads in the region
                let intervals = records
                    .into_iter()
                    .filter(|r| r.is_ok())
                    .map(|r| r.unwrap())
                    .map(|r| {
                        IntervalMaker::new(
                            r,
                            &header,
                            &chromsizes_refid,
                            &self.filter,
                            self.use_fragment,
                            None,
                        )
                    })
                    .map(|i| i.coords())
                    .filter(|c| c.is_some())
                    .map(|c| c.unwrap())
                    .map(|i| Iv {
                        start: i.0,
                        stop: i.1,
                        val: 1,
                    })
                    .collect::<Vec<Iv>>();

                // Create a lapper from the intervals
                let mut lapper = Lapper::new(intervals);

                // Iterate over bins within the given region and count the reads using the lapper
                let mut bin_counts: Vec<Iv> = Vec::new();

                // Generate the bins for the region (bin_size)
                let region_interval = region.interval();
                let region_start = region_interval
                    .start()
                    .expect("Error getting interval start")
                    .get();
                let region_end = region_interval
                    .end()
                    .expect("Error getting interval end")
                    .get();
                let bin_starts = (region_start..region_end).step_by(self.bin_size as usize); // Corrected to include the end in the range
                let bin_ends = bin_starts
                    .clone()
                    .skip(1);
                let n_bins =
                    ((region_end - region_start) + self.bin_size as usize) / self.bin_size as usize; // Adjusted to correctly calculate the number of bins when the range is not a perfect multiple of bin_size

                // Iterate over the bins
                for (start, end) in bin_starts.zip(bin_ends) {
                    let count = lapper.count(start as usize, end as usize);

                    // Ensure that the score is not the same as the previous entry. Need to just collapse these to save memory
                    match bin_counts.last_mut() {
                        Some(entry) => {
                            if entry.val == count as u32 {
                                entry.stop = end;
                                continue;
                            } else {
                                bin_counts.push(Iv {
                                    start: start,
                                    stop: end,
                                    val: count as u32,
                                });
                                continue;
                            }

                        }
                        None => {
                            bin_counts.push(Iv {
                                start: start,
                                stop: end,
                                val: count as u32,
                            });
                            continue;
                        }
                    };
                }

                (region.name().to_owned().to_string(), bin_counts)
            })
            .fold(
                || HashMap::new(),
                |mut acc, (chrom, intervals)| {
                    acc.entry(chrom).or_insert_with(Vec::new).extend(intervals);
                    acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut acc, mut map| {
                    for (key, value) in map.drain() {
                        acc.entry(key).or_insert_with(Vec::new).extend(value);
                    }
                    acc
                },
            );

        info!("Pileup complete");
        info!("Read filtering statistics: {}", self.filter.stats());

        // Assuming here that the pileups are returned in the correct order
        Ok(pileup)
    }

    fn pileup_to_polars(&self) -> Result<DataFrame> {
        // Get the pileup
        let pileup = self.pileup().expect("Error getting pileup");

        // Create a DataFrame from the pileup
        let (chrom, start, end, score) = pileup
            .into_par_iter()
            .map(|(chrom, intervals)| {
                let chrom = vec![chrom; intervals.len()];
                let start = intervals
                    .iter()
                    .map(|i| i.start as u64)
                    .collect::<Vec<u64>>();
                let end = intervals
                    .iter()
                    .map(|i| i.stop as u64)
                    .collect::<Vec<u64>>();
                let score = intervals.iter().map(|i| i.val as u32).collect::<Vec<u32>>();

                (chrom, start, end, score)
            })
            .reduce(
                || (vec![], vec![], vec![], vec![]),
                |(mut chrom, mut start, mut end, mut score), (chrom_, start_, end_, score_)| {
                    chrom.extend(chrom_);
                    start.extend(start_);
                    end.extend(end_);
                    score.extend(score_);

                    (chrom, start, end, score)
                },
            );

        Ok(DataFrame::new(vec![
            Series::new("chrom", chrom),
            Series::new("start", start),
            Series::new("end", end),
            Series::new("score", score),
        ])?)
    }

    fn pileup_normalised(&self) -> Result<DataFrame> {
        
        let mut df = self.pileup_to_polars().expect("Error getting pileup");
        let n_total_counts = df.column("score")?.sum::<u32>().unwrap() as u64;
        info!("Total counts: {}", n_total_counts);

        // Take care of the normalization - this will only work for CPM correctly at the moment
        let norm_factor =
            self.norm_method
                .scale_factor(self.scale_factor, self.bin_size, n_total_counts);

        // Normalise the score column
        let score = df.column("score")?.u32()? * norm_factor;
        df.replace("score", score)?;

        Ok(df)
    }

    pub fn to_bedgraph(&self, outfile: PathBuf) -> Result<()> {

        info!("Writing bedgraph file to {}", outfile.display());
        let mut file =
            std::fs::File::create(outfile).expect("could not create output bedgraph file");

        let mut df = self.pileup_normalised().expect("Error getting pileup");
        df.sort(["chrom", "start"], Default::default())?;

        CsvWriter::new(&mut file)
            .include_header(false)
            .with_separator(b'\t')
            .finish(&mut df)?;

        Ok(())
    }
    pub fn to_bigwig(&self, outfile: PathBuf) -> Result<()> {
        // Write chromosome sizes to file
        let bam_details = BamStats::new(self.file_path.clone())?;
        let chromsizes_file = tempfile::NamedTempFile::new()?;
        let chromsizes_path = chromsizes_file.path();
        bam_details.write_chromsizes(chromsizes_path.to_path_buf())?;

        // Write the bedgraph file
        let bedgraph_file = tempfile::NamedTempFile::new()?;
        let bedgraph_path = bedgraph_file.path();
        self.to_bedgraph(bedgraph_path.to_path_buf())?;

        // Convert the bedgraph file to bigwig
        info!("Converting bedGraph to BigWig file");
        let bw_command_output = std::process::Command::new("bedGraphToBigWig")
            .arg(bedgraph_path)
            .arg(chromsizes_path)
            .arg(outfile.clone())
            .output()
            .expect("Failed to convert bedgraph to bigwig");

        // Check if the command was successful
        if !bw_command_output.status.success() {
            log::error!("Failed to convert bedgraph to bigwig");
            log::error!("{}", String::from_utf8_lossy(&bw_command_output.stderr));
            std::process::exit(1);
        } else {
            info!("BigWig file written to {}", outfile.display());
        }

        Ok(())
    }
}
