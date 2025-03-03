//! This module implements pileup generation from one or multiple BAM files
//! and converts the resulting signal into bedGraph and BigWig formats. The
//! implementation is parallelized (using Rayon) and uses several libraries to
//! process genomic intervals and to normalize and aggregate counts.
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::ops::Deref;
use std::path::PathBuf;
use std::process::Command;

use ahash::HashMap;
use anyhow::{Context, Result};
use indicatif::ParallelProgressIterator;
use log::{debug, error, info};
use noodles::{bam, core};
use polars::lazy::dsl::{col, cols, mean_horizontal, sum_horizontal};
use polars::prelude::*;
use rayon::prelude::*;
use rust_lapper::Lapper;
use tempfile;

use crate::filter::BamReadFilter;
use crate::intervals::IntervalMaker;
use crate::utils::{get_bam_header, progress_bar, BamStats, Iv, NormalizationMethod};

/// Write a DataFrame as a bedGraph file (tab-separated, no header).
fn write_bedgraph(mut df: DataFrame, outfile: PathBuf) -> Result<()> {
    // Sort by chromosome and start position.
    df = df.sort(["chrom", "start"], Default::default())?;
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
fn collapse_equal_bins(
    df:  DataFrame,
    score_columns: Option<Vec<PlSmallStr>>,
) -> Result<DataFrame> {
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

/// Convert a bedGraph file to BigWig format using the external command
/// `bedGraphToBigWig`.
fn convert_bedgraph_to_bigwig(
    bedgraph_path: &std::path::Path,
    chromsizes_path: &std::path::Path,
    outfile: &PathBuf,
) -> Result<()> {
    let output = Command::new("bedGraphToBigWig")
        .arg(bedgraph_path)
        .arg(chromsizes_path)
        .arg(outfile)
        .output()
        .context("Failed to execute bedGraphToBigWig")?;

    if !output.status.success() {
        error!("Error converting bedGraph to BigWig:");
        error!("{}", String::from_utf8_lossy(&output.stderr));
        anyhow::bail!("Conversion to BigWig failed");
    } else {
        info!("BigWig file successfully written to {}", outfile.display());
    }
    Ok(())
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
    pub fn new(
        file_path: PathBuf,
        bin_size: u64,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        use_fragment: bool,
        filter: BamReadFilter,
        colapse: bool,
    ) -> Self {
        Self {
            file_path,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
            filter,
            collapse: colapse,
        }
    }

    /// Generate pileup intervals for the BAM file.
    ///
    /// Returns a map from chromosome names to a vector of intervals (with
    /// start, stop and count) for each genomic chunk.
    fn pileup(&self) -> Result<HashMap<String, Vec<Iv>>> {
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let genomic_chunks = bam_stats.genome_chunks(self.bin_size)?;
        let chromsizes_refid = bam_stats.chromsizes_ref_id()?;
        let n_total_chunks = genomic_chunks.len();

        info!("{}", self);
        info!("Processing {} genomic chunks", n_total_chunks);

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
                let header = get_bam_header(self.file_path.clone())?;

                // Query for reads overlapping the region.
                let records = reader.query(&header, &region)?;

                // Create intervals from each read that passes filtering.
                let intervals: Vec<Iv> = records
                    .into_iter()
                    .filter_map(|record| record.ok())
                    .filter_map(|record| {
                        IntervalMaker::new(
                            record,
                            &header,
                            &chromsizes_refid,
                            &self.filter,
                            self.use_fragment,
                            None,
                        )
                        .coords()
                        .map(|(s, e)| Iv {
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
                    let count = lapper.count(start as usize, end as usize);

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
                || HashMap::default(),
                |mut acc, result: Result<(String, Vec<Iv>)>| {
                    if let Ok((chrom, intervals)) = result {
                        acc.entry(chrom).or_insert_with(Vec::new).extend(intervals);
                    }
                    acc
                },
            )
            .reduce(
                || HashMap::default(),
                |mut acc, map| {
                    for (key, mut value) in map {
                        acc.entry(key).or_insert_with(Vec::new).append(&mut value);
                    }
                    acc
                },
            );

        info!("Pileup complete");
        info!("Read filtering statistics: {}", self.filter.stats());

        Ok(pileup)
    }

    /// Convert the pileup intervals into a Polars DataFrame.
    ///
    /// The DataFrame will have columns "chrom", "start", "end" and "score".
    fn pileup_to_polars(&self) -> Result<DataFrame> {
        let pileup = self.pileup()?;

        // Process each chromosome in parallel and combine into column vectors.
        let (chroms, starts, ends, scores) = pileup
            .into_par_iter()
            .map(|(chrom, intervals)| {
                let chrom_vec = vec![chrom; intervals.len()];
                let start_vec = intervals
                    .iter()
                    .map(|iv| iv.start as u64)
                    .collect::<Vec<_>>();
                let end_vec = intervals
                    .iter()
                    .map(|iv| iv.stop as u64)
                    .collect::<Vec<_>>();
                let score_vec = intervals.iter().map(|iv| iv.val as u32).collect::<Vec<_>>();
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

        // Build the DataFrame.
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
        // Get the total counts across all bins.
        let n_total_counts: u64 = df
            .column("score")?
            .sum_reduce()?
            .as_any_value()
            .try_extract()?;
        info!("Total counts: {}", n_total_counts);

        // Compute the normalization factor.
        let norm_factor =
            self.norm_method
                .scale_factor(self.scale_factor, self.bin_size, n_total_counts);
        // Multiply the score column by the normalization factor.
        let norm_scores = df.column("score")?.u32()? * norm_factor;
        df.replace("score", norm_scores)?;
        Ok(df)
    }

    /// Write the normalized pileup as a bedGraph file.
    pub fn to_bedgraph(&self, outfile: PathBuf) -> Result<()> {
        info!("Writing bedGraph file to {}", outfile.display());
        let df = self.pileup_normalised()?;
        write_bedgraph(df, outfile)
    }

    /// Write the normalized pileup as a BigWig file.
    ///
    /// This function writes a temporary bedGraph file and converts it to BigWig.
    pub fn to_bigwig(&self, outfile: PathBuf) -> Result<()> {
        let bam_stats = BamStats::new(self.file_path.clone())?;
        let chromsizes_file = tempfile::NamedTempFile::new()?;
        let chromsizes_path = chromsizes_file.path();
        bam_stats.write_chromsizes(chromsizes_path.to_path_buf())?;

        let bedgraph_file = tempfile::NamedTempFile::new()?;
        let bedgraph_path = bedgraph_file.path();
        self.to_bedgraph(bedgraph_path.to_path_buf())?;

        info!("Converting bedGraph to BigWig file");
        convert_bedgraph_to_bigwig(bedgraph_path, chromsizes_path, &outfile)
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
        if collapse.iter().all(|x| *x == false) {
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
        let df = df.join(
            scores_df,
            [col("chrom"), col("start"), col("end")],
            [col("chrom"), col("start"), col("end")],
            JoinArgs::new(JoinType::Full),
        ).collect()?;
        let df = collapse_equal_bins(df, Some(names))?;

        Ok(df)
    }
}
