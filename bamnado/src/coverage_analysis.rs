//! # Coverage Analysis Module
//!
//! This module implements the core logic for generating coverage tracks (pileups) from
//! one or multiple BAM files. It supports:
//! *   Parallel processing of BAM files using Rayon.
//! *   Generation of bedGraph and BigWig output formats.
//! *   Normalization of coverage signals (RPKM, CPM, etc.).
//! *   Handling of genomic intervals and binning.

use std::fmt::{Display, Formatter, Result as FmtResult};
use std::path::PathBuf;

use crate::bam_utils::{BamStats, Iv, get_bam_header, progress_bar};
use crate::genomic_intervals::{IntervalMaker, Shift, Truncate};
use crate::read_filter::BamReadFilter;
use crate::signal_normalization::NormalizationMethod;
use ahash::HashMap;
use anyhow::{Context, Result};
#[cfg(test)]
use bigtools::BigWigRead;
use bigtools::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};
use indicatif::ParallelProgressIterator;
use log::{debug, error, info};
use noodles::bam;
use polars::lazy::dsl::col;
use polars::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use rust_lapper::Lapper;
use tempfile;

/// Write a DataFrame as a bedGraph file (tab-separated, no header).
fn write_bedgraph(
    df: DataFrame,
    outfile: PathBuf,
    ignore_scaffold_chromosomes: bool,
) -> Result<()> {
    // If ignore_scaffold_chromosomes is true, filter out scaffold chromosomes
    let mut df = match ignore_scaffold_chromosomes {
        true => {
            let pattern = Regex::new("(chrUn.*|.*alt|.*random)").unwrap();
            let mask = df
                .column("chrom")?
                .str()?
                .iter()
                .flatten()
                .map(|s| !pattern.is_match(s))
                .collect::<BooleanChunked>();
            df.filter(&mask)?
        }
        false => df,
    };

    df.sort_in_place(["chrom", "start"], SortMultipleOptions::default())?;

    // Select only the required columns in the correct order for bedGraph format
    let mut df = df.select(["chrom", "start", "end", "score"])?;

    let mut file = std::fs::File::create(outfile)?;
    CsvWriter::new(&mut file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(&mut df)?;

    Ok(())
}

/// Collapse adjacent bins that have equal chromosome and score values.
///
/// This function reduces the number of rows by merging adjacent bins with the
/// same score, which can make the output bedGraph file smaller.
fn collapse_equal_bins(df: DataFrame, score_columns: Option<Vec<PlSmallStr>>) -> Result<DataFrame> {
    let mut df = df.sort(["chrom", "start"], Default::default())?;
    let shifted_chromosome = df.column("chrom")?.shift(1);
    let same_chrom = df.column("chrom")?.equal(&shifted_chromosome)?;

    // For each row check that all the scores are the same in the shifted score columns as the original score columns
    let score_columns = score_columns.unwrap_or_else(|| vec!["score".into()]);

    let scores = df.select(score_columns.clone())?;
    let scores_shifted = df.select(score_columns.clone())?.shift(1);
    let mut same_scores = Vec::new();
    for (score, score_shifted) in scores.iter().zip(scores_shifted.iter()) {
        let same_score = score.equal(score_shifted)?;
        same_scores.push(same_score);
    }

    // Sum the boolean columns to get a single boolean column
    let same_chrom_and_score = same_chrom
        & same_scores
            .iter()
            .fold(same_scores[0].clone(), |acc, x| acc & x.clone());

    // Compute a cumulative sum to define groups of identical rows.
    let group = cum_sum(&same_chrom_and_score.into_series(), false)?
        .with_name("groups".into())
        .into_column();

    let df = df
        .with_column(group)?
        .clone()
        .lazy()
        .group_by(["groups"])
        .agg(&[
            col("chrom").first(),
            col("start").min(),
            col("end").max(),
            col("score").sum(),
        ])
        .collect()?;

    Ok(df)
}

pub struct BamPileup {
    // Path to the BAM file.
    file_path: PathBuf,
    // Bin size for counting reads.
    bin_size: u64,
    // Normalization method for the pileup signal.
    norm_method: NormalizationMethod,
    // Scaling factor for the normalization method.
    scale_factor: f32,
    // Use read fragments for counting instead of individual reads.
    use_fragment: bool,
    // Filter for reads to include in the pileup.
    filter: BamReadFilter,
    // Merge adjacent bins with equal scores.
    collapse: bool,
    // Ignore scaffold chromosomes in the pileup.
    ignore_scaffold_chromosomes: bool,
    // Optional shift to apply to the pileup intervals.
    shift: Option<Shift>,
    // Optional truncation to apply to the pileup intervals.
    truncate: Option<Truncate>,
}

impl Display for BamPileup {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        writeln!(f, "Pileup Settings:")?;
        writeln!(f, "BAM file: {}", self.file_path.display())?;
        writeln!(f, "Bin size: {}", self.bin_size)?;
        writeln!(f, "Normalization method: {:?}", self.norm_method)?;
        writeln!(f, "Scaling factor: {}", self.scale_factor)?;
        writeln!(f, "Using fragment for counting: {}", self.use_fragment)?;
        write!(f, "Filtering using: \n{}\n", self.filter)
    }
}

impl BamPileup {
    /// Create a new [`BamPileup`].
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        file_path: PathBuf,
        bin_size: u64,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        use_fragment: bool,
        filter: BamReadFilter,
        collapse_intervals: bool,
        ignore_scaffold_chromosomes: bool,
        shift: Option<Shift>,
        truncate: Option<Truncate>,
    ) -> Self {
        Self {
            file_path,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
            filter,
            collapse: collapse_intervals,
            ignore_scaffold_chromosomes,
            shift,
            truncate,
        }
    }

    pub fn normalization_factor(&self) -> f64 {
        // Calculate the normalization factor based on the normalization method.
        let filter_stats = self.filter.stats();
        let n_reads = filter_stats.n_reads_after_filtering();
        let normalization_factor =
            self.norm_method
                .scale_factor(self.scale_factor, self.bin_size, n_reads);
        // Log the normalization factor.
        info!("Normalization factor: {normalization_factor}");
        normalization_factor
    }

    /// Generate pileup for a single chromosome.
    /// Returns a vector of intervals (with start, stop and count) for the chromosome.
    pub fn pileup_chromosome(&self, chrom: &str) -> Result<Vec<Iv>> {
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let genomic_chunks = bam_stats.chromosome_chunks(chrom, self.bin_size)?;
        let chromsizes_refid = bam_stats.chromsizes_ref_id()?;
        let n_total_chunks = genomic_chunks.len();

        let header = get_bam_header(self.file_path.clone())?;
        // Process each genomic chunk in parallel.
        let pileup: Vec<Iv> = genomic_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                n_total_chunks as u64,
                format!("Performing pileup for chromosome {chrom}"),
            ))
            .map(|region| {
                // Each thread creates its own BAM reader.
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(self.file_path.clone())
                    .context("Failed to open BAM file")?;
                // Query for reads overlapping the region.
                let records = reader.query(&header, &region)?;
                // Create intervals from each read that passes filtering.
                let intervals: Vec<Iv> = records
                    .records()
                    .filter_map(|record| record.ok())
                    .filter_map(|record| {
                        IntervalMaker::new(
                            record,
                            &header,
                            &chromsizes_refid,
                            &self.filter,
                            self.use_fragment,
                            self.shift,
                            self.truncate,
                        )
                        .coords()
                        .map(|(s, e, _dtlen)| Iv {
                            start: s,
                            stop: e,
                            val: 1,
                        })
                    })
                    .collect();
                // Use a Lapper to count overlapping intervals in bins.
                let lapper = Lapper::new(intervals);
                let region_interval = region.interval();
                let region_start = region_interval
                    .start()
                    .context("Failed to get region start")?
                    .get();
                let region_end = region_interval
                    .end()
                    .context("Failed to get region end")?
                    .get();
                let mut bin_counts: Vec<Iv> = Vec::new();
                let mut start = region_start;
                while start < region_end {
                    let end = (start + self.bin_size as usize).min(region_end);
                    let count = lapper.count(start, end);
                    match self.collapse {
                        true => {
                            if let Some(last) = bin_counts.last_mut() {
                                if last.val == count as u32 {
                                    last.stop = end;
                                } else {
                                    bin_counts.push(Iv {
                                        start,
                                        stop: end,
                                        val: count as u32,
                                    });
                                }
                            } else {
                                bin_counts.push(Iv {
                                    start,
                                    stop: end,
                                    val: count as u32,
                                });
                            }
                        }
                        false => {
                            bin_counts.push(Iv {
                                start,
                                stop: end,
                                val: count as u32,
                            });
                        }
                    }
                    start = end;
                }
                Ok(bin_counts)
            })
            // Combine the results from parallel threads.
            .fold(Vec::new, |mut acc, result: Result<Vec<Iv>>| {
                if let Ok(mut intervals) = result {
                    acc.append(&mut intervals);
                }
                acc
            })
            .reduce(Vec::new, |mut acc, mut vec| {
                acc.append(&mut vec);
                acc
            });
        Ok(pileup)
    }

    /// Generate pileup intervals for the BAM file.
    ///
    /// Returns a map from chromosome names to a vector of intervals (with
    /// start, stop and count) for each genomic chunk.
    fn pileup_genome(&self) -> Result<HashMap<String, Vec<Iv>>> {
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let genomic_chunks = bam_stats.genome_chunks(self.bin_size)?;
        let chromsizes_refid = bam_stats.chromsizes_ref_id()?;
        let n_total_chunks = genomic_chunks.len();

        info!("{self}");
        info!("Processing {n_total_chunks} genomic chunks");

        let header = get_bam_header(self.file_path.clone())?;

        // Process each genomic chunk in parallel.
        let pileup = genomic_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                n_total_chunks as u64,
                "Performing pileup".to_string(),
            ))
            .map(|region| {
                // Each thread creates its own BAM reader.
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(self.file_path.clone())
                    .context("Failed to open BAM file")?;

                // Query for reads overlapping the region.
                let records = reader.query(&header, &region)?;

                // Create intervals from each read that passes filtering.
                let intervals: Vec<Iv> = records
                    .records()
                    .filter_map(|record| record.ok())
                    .filter_map(|record| {
                        IntervalMaker::new(
                            record,
                            &header,
                            &chromsizes_refid,
                            &self.filter,
                            self.use_fragment,
                            self.shift,
                            self.truncate,
                        )
                        .coords()
                        .map(|(s, e, _dtlen)| Iv {
                            start: s,
                            stop: e,
                            val: 1,
                        })
                    })
                    .collect();

                // Use a Lapper to count overlapping intervals in bins.
                let lapper = Lapper::new(intervals);
                let region_interval = region.interval();
                let region_start = region_interval
                    .start()
                    .context("Failed to get region start")?
                    .get();
                let region_end = region_interval
                    .end()
                    .context("Failed to get region end")?
                    .get();

                let mut bin_counts: Vec<rust_lapper::Interval<usize, u32>> = Vec::new();
                let mut start = region_start;
                while start < region_end {
                    let end = (start + self.bin_size as usize).min(region_end);
                    let count = lapper.count(start, end);

                    match self.collapse {
                        true => {
                            if let Some(last) = bin_counts.last_mut() {
                                if last.val == count as u32 {
                                    last.stop = end;
                                } else {
                                    bin_counts.push(rust_lapper::Interval {
                                        start,
                                        stop: end,
                                        val: count as u32,
                                    });
                                }
                            } else {
                                bin_counts.push(rust_lapper::Interval {
                                    start,
                                    stop: end,
                                    val: count as u32,
                                });
                            }
                        }
                        false => {
                            bin_counts.push(rust_lapper::Interval {
                                start,
                                stop: end,
                                val: count as u32,
                            });
                        }
                    }
                    start = end;
                }

                Ok((region.name().to_owned().to_string(), bin_counts))
            })
            // Combine the results from parallel threads.
            .fold(
                HashMap::<String, Vec<Iv>>::default,
                |mut acc, result: Result<(String, Vec<Iv>)>| {
                    if let Ok((chrom, intervals)) = result {
                        acc.entry(chrom).or_default().extend(intervals);
                    }
                    acc
                },
            )
            .reduce(HashMap::default, |mut acc, map| {
                for (key, mut value) in map {
                    acc.entry(key).or_default().append(&mut value);
                }
                acc
            });

        info!("Pileup complete");
        info!("Read filtering statistics: {}", self.filter.stats());

        Ok(pileup)
    }

    /// Convert the pileup intervals into a Polars DataFrame.
    ///
    /// The DataFrame will have columns "chrom", "start", "end" and "score".
    fn pileup_to_polars(&self) -> Result<DataFrame> {
        let pileup = self.pileup_genome()?;

        // Process each chromosome in parallel and combine into column vectors.
        let (chroms, starts, ends, scores) = pileup
            .into_par_iter()
            .map(|(chrom, intervals)| {
                let n = intervals.len();
                // Preallocate with exact capacity.
                let mut start_vec = Vec::with_capacity(n);
                let mut end_vec = Vec::with_capacity(n);
                let mut score_vec = Vec::with_capacity(n);
                // Process each interval in one pass.
                for iv in intervals {
                    start_vec.push(iv.start as u64);
                    end_vec.push(iv.stop as u64);
                    score_vec.push(iv.val);
                }
                // Create the chrom vector by repeating the chromosome name.
                let chrom_vec = vec![chrom; n];
                (chrom_vec, start_vec, end_vec, score_vec)
            })
            .reduce(
                || (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                |(mut chrom_a, mut start_a, mut end_a, mut score_a),
                 (chrom_b, start_b, end_b, score_b)| {
                    chrom_a.extend(chrom_b);
                    start_a.extend(start_b);
                    end_a.extend(end_b);
                    score_a.extend(score_b);
                    (chrom_a, start_a, end_a, score_a)
                },
            );

        // Build the DataFrame using the combined vectors.
        let df = DataFrame::new(vec![
            Column::new("chrom".into(), chroms),
            Column::new("start".into(), starts),
            Column::new("end".into(), ends),
            Column::new("score".into(), scores),
        ])?;
        Ok(df)
    }

    /// Normalize the pileup signal using the provided normalization method.
    fn pileup_normalised(&self) -> Result<DataFrame> {
        let mut df = self.pileup_to_polars()?;
        let norm_scores =
            df.column("score")?.cast(&DataType::Float64)? * self.normalization_factor();
        df.with_column(norm_scores)?;
        Ok(df)
    }

    /// Write the normalized pileup as a bedGraph file.
    ///
    /// # Arguments
    ///
    /// * `outfile` - The path to the output bedGraph file.
    pub fn to_bedgraph(&self, outfile: PathBuf) -> Result<()> {
        info!("Writing bedGraph file to {}", outfile.display());
        let df = self.pileup_normalised()?;
        write_bedgraph(df, outfile, self.ignore_scaffold_chromosomes)?;
        Ok(())
    }

    /// Write the normalized pileup as a BigWig file.
    ///
    /// This function writes a temporary bedGraph file and converts it to BigWig.
    ///
    /// # Arguments
    ///
    /// * `outfile` - The path to the output BigWig file.
    pub fn to_bigwig(&self, outfile: PathBuf) -> Result<()> {
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let chromsizes_file = tempfile::NamedTempFile::new()?;
        let chromsizes_path = chromsizes_file.path();
        bam_stats.write_chromsizes(chromsizes_path.to_path_buf())?;

        let bedgraph_file = tempfile::NamedTempFile::new()?;
        let bedgraph_path = bedgraph_file.path();
        self.to_bedgraph(bedgraph_path.to_path_buf())?;

        let args = bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs {
            bedgraph: bedgraph_path.to_path_buf().to_string_lossy().to_string(),
            chromsizes: chromsizes_path.to_path_buf().to_string_lossy().to_string(),
            output: outfile.to_string_lossy().to_string(),
            parallel: "auto".to_string(),
            single_pass: true,
            write_args: bigtools::utils::cli::BBIWriteArgs {
                nthreads: 6,
                nzooms: 10,
                uncompressed: false,
                sorted: "all".to_string(),
                zooms: None,
                block_size: DEFAULT_BLOCK_SIZE,
                items_per_slot: DEFAULT_ITEMS_PER_SLOT,
                inmemory: false,
            },
        };

        let result = bigtools::utils::cli::bedgraphtobigwig::bedgraphtobigwig(args);
        if let Err(e) = result {
            error!("Error converting bedGraph to BigWig: {e}");
            anyhow::bail!("Conversion to BigWig failed");
        }
        // Clean up temporary files.
        std::fs::remove_file(bedgraph_path)?;
        std::fs::remove_file(chromsizes_path)?;
        // Log success.
        info!("BigWig file successfully written to {}", outfile.display());
        Ok(())
    }
}

/// A struct for processing multiple BAM files in parallel.
/// This struct is useful for combining pileups from multiple samples.
/// The pileups are processed in parallel and combined into a single DataFrame.
/// The resulting DataFrame can be written to a bedGraph or BigWig file.
pub struct MultiBamPileup {
    // BamPileups to process in parallel.
    bam_pileups: Vec<BamPileup>,
}

impl MultiBamPileup {
    /// Create a new [`MultiBamPileup`] from a vector of [`BamPileup`]s.
    ///
    /// # Arguments
    ///
    /// * `pileups` - A vector of `BamPileup` structs.
    pub fn new(pileups: Vec<BamPileup>) -> Self {
        Self {
            bam_pileups: pileups,
        }
    }

    fn _validate_bam_pileups(&self) -> Result<()> {
        // Check that all BAM pileups have the same bin size.
        let bin_sizes: Vec<u64> = self.bam_pileups.iter().map(|p| p.bin_size).collect();
        if bin_sizes.iter().all(|&x| x == bin_sizes[0]) {
            debug!("All BAM pileups have the same bin size: {}", bin_sizes[0]);
        } else {
            anyhow::bail!("All BAM pileups must have the same bin size");
        }

        // Check that all BAM pileups have the collapse option set to false
        let collapse: Vec<bool> = self.bam_pileups.iter().map(|p| p.collapse).collect();
        if collapse.iter().all(|x| !(*x)) {
            debug!("All BAM pileups have collapse set to false");
        } else {
            anyhow::bail!("All BAM pileups must have collapse set to false");
        }

        Ok(())
    }

    /// Generate pileup intervals for all BAM files.
    ///
    /// Returns a map from chromosome names to a vector of intervals (with
    /// start, stop and count) for each genomic chunk.
    fn pileup(&self) -> Result<DataFrame> {
        self._validate_bam_pileups()?;

        let mut dfs = Vec::new();
        let names: Vec<PlSmallStr> = self
            .bam_pileups
            .iter()
            .map(|p| p.file_path.clone().to_str().unwrap().into())
            .collect::<Vec<_>>();

        for pileup in &self.bam_pileups {
            let df = pileup.pileup_to_polars()?;
            dfs.push(df);
        }

        let df = dfs[0]
            .clone()
            .lazy()
            .rename(["score"], [names[0].clone()], true);
        let scores = dfs
            .iter()
            .skip(1)
            .zip(names.iter().skip(1))
            .map(|(df, n)| df.column("score").unwrap().clone().with_name(n.clone()))
            .collect::<Vec<Column>>();
        let scores_df = DataFrame::new(scores)?.lazy();
        let df = df
            .join(
                scores_df,
                [col("chrom"), col("start"), col("end")],
                [col("chrom"), col("start"), col("end")],
                JoinArgs::new(JoinType::Full),
            )
            .collect()?;
        let df = collapse_equal_bins(df, Some(names))?;

        Ok(df)
    }

    /// Writes the multi-BAM pileup to a TSV file.
    ///
    /// # Arguments
    ///
    /// * `outfile` - The path to the output TSV file.
    pub fn to_tsv(&self, outfile: &PathBuf) -> Result<()> {
        let mut df = self.pileup()?;
        let mut file = std::fs::File::create(outfile)?;
        CsvWriter::new(&mut file)
            .include_header(true)
            .with_separator(b'\t')
            .finish(&mut df)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_filter::BamReadFilter;
    use crate::signal_normalization::NormalizationMethod;
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Helper function to create a test BamPileup instance
    fn create_test_bam_pileup() -> BamPileup {
        BamPileup::new(
            PathBuf::from("test.bam"),
            1000,
            NormalizationMethod::Raw,
            1.0,
            false,
            BamReadFilter::default(),
            false,
            false,
            None,
            None,
        )
    }

    fn get_test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("Failed to get workspace root")
            .join("test/data")
    }

    #[test]
    fn test_bam_pileup_new() {
        let file_path = PathBuf::from("test.bam");
        let pileup = BamPileup::new(
            file_path.clone(),
            1000,
            NormalizationMethod::RPKM,
            2.0,
            true,
            BamReadFilter::default(),
            true,
            false,
            None,
            None,
        );

        assert_eq!(pileup.file_path, file_path);
        assert_eq!(pileup.bin_size, 1000);
        assert_eq!(pileup.scale_factor, 2.0);
        assert!(pileup.use_fragment);
        assert!(pileup.collapse);
    }

    #[test]
    fn test_write_bedgraph() {
        let temp_dir = TempDir::new().unwrap();
        let output_path = temp_dir.path().join("test.bedgraph");

        // Create test DataFrame
        let df = DataFrame::new(vec![
            Column::new("chrom".into(), vec!["chr1", "chr1", "chr2"]),
            Column::new("start".into(), vec![0u64, 1000u64, 0u64]),
            Column::new("end".into(), vec![1000u64, 2000u64, 1000u64]),
            Column::new("score".into(), vec![10u32, 20u32, 15u32]),
        ])
        .unwrap();

        let result = write_bedgraph(df, output_path.clone(), true);
        assert!(result.is_ok());
        assert!(output_path.exists());

        // Verify the file was written correctly
        let contents = std::fs::read_to_string(output_path).unwrap();
        let lines: Vec<&str> = contents.trim().split('\n').collect();
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn test_collapse_equal_bins() {
        // Create DataFrame with adjacent bins having same scores
        let df = DataFrame::new(vec![
            Column::new("chrom".into(), vec!["chr1", "chr1", "chr1", "chr1"]),
            Column::new("start".into(), vec![0u64, 1000u64, 2000u64, 3000u64]),
            Column::new("end".into(), vec![1000u64, 2000u64, 3000u64, 4000u64]),
            Column::new("score".into(), vec![10u32, 10u32, 20u32, 20u32]),
        ])
        .unwrap();

        let result = collapse_equal_bins(df, None).unwrap();

        // Should collapse adjacent bins with same scores
        assert!(result.height() < 4);

        // Verify the collapsed bins have correct aggregated values
        let score_col = result.column("score").unwrap();
        let score = score_col.sum_reduce().expect("Failed to sum scores");
        assert_eq!(score.as_any_value().try_extract::<u32>().unwrap(), 60);
    }

    #[test]
    fn test_bam_pileup_display() {
        let pileup = create_test_bam_pileup();
        let display_str = format!("{pileup}");

        assert!(display_str.contains("Pileup Settings:"));
        assert!(display_str.contains("BAM file: test.bam"));
        assert!(display_str.contains("Bin size: 1000"));
        assert!(display_str.contains("Normalization method:"));
        assert!(display_str.contains("Scaling factor: 1"));
        assert!(display_str.contains("Using fragment for counting: false"));
    }

    #[test]
    fn test_bam_to_bigwig_raw_scale() {
        let test_data_dir = get_test_data_dir();
        let temp_dir = TempDir::new().unwrap();
        let bam_file = test_data_dir.join("test.bam");
        let output_path = temp_dir.path().join("test.bw");
        let pileup = BamPileup::new(
            bam_file,
            1000,
            NormalizationMethod::Raw,
            1.0,
            false,
            BamReadFilter::default(),
            true,
            false,
            None,
            None,
        );
        let result = pileup.to_bigwig(output_path.clone());
        result.expect("Failed to create BigWig file");

        assert!(output_path.exists(), "BigWig file was not created");
        info!("BigWig file created at: {}", output_path.display());

        let bw = BigWigRead::open_file(output_path.clone());
        assert!(bw.is_ok(), "Failed to open BigWig file");
        let mut bw = bw.unwrap();

        let stats = bw.get_summary().expect("Failed to get BigWig summary");

        assert_eq!(stats.min_val, 0.0);
        assert_eq!(stats.max_val, 117.0);
    }

    #[test]
    fn test_bam_to_bigwig_cpm_scale() {
        let test_data_dir = get_test_data_dir();
        let temp_dir = TempDir::new().unwrap();
        let bam_file = test_data_dir.join("test.bam");
        let output_path = temp_dir.path().join("test_cpm.bw");
        let pileup = BamPileup::new(
            bam_file,
            1000,
            NormalizationMethod::CPM,
            1.0,
            false,
            BamReadFilter::default(),
            true,
            false,
            None,
            None,
        );
        let result = pileup.to_bigwig(output_path.clone());
        result.expect("Failed to create BigWig file");

        assert!(output_path.exists(), "BigWig file was not created");
        info!("BigWig file created at: {}", output_path.display());

        let bw = BigWigRead::open_file(output_path.clone());
        assert!(bw.is_ok(), "Failed to open BigWig file");
        let mut bw = bw.unwrap();

        let stats = bw.get_summary().expect("Failed to get BigWig summary");

        println!("BigWig stats: {stats:?}");
        assert_eq!(stats.min_val, 0.0);
        assert_eq!(stats.max_val, 10134.2568359375);
    }

    #[test]
    fn test_bam_to_bigwig_rpkm_scale() {
        let test_data_dir = get_test_data_dir();
        let temp_dir = TempDir::new().unwrap();
        let bam_file = test_data_dir.join("test.bam");
        let output_path = temp_dir.path().join("test_rpkm.bw");
        let pileup = BamPileup::new(
            bam_file,
            1000,
            NormalizationMethod::RPKM,
            1.0,
            false,
            BamReadFilter::default(),
            false,
            false,
            None,
            None,
        );
        let result = pileup.to_bigwig(output_path.clone());
        result.expect("Failed to create BigWig file");

        assert!(output_path.exists(), "BigWig file was not created");
        info!("BigWig file created at: {}", output_path.display());

        let bw = BigWigRead::open_file(output_path.clone());
        assert!(bw.is_ok(), "Failed to open BigWig file");
        let mut bw = bw.unwrap();

        let stats = bw.get_summary().expect("Failed to get BigWig summary");

        println!("BigWig stats: {stats:?}");
        assert_eq!(stats.min_val, 0.0);
        assert_eq!(stats.max_val, 10134.2568359375);
    }
}
