//! # Count Module
//!
//! This module implements pileup generation from one or multiple BAM files
//! and converts the resulting signal into bedGraph and BigWig formats.
//!
//! It supports:
//! *   Parallel processing using Rayon.
//! *   Normalization of counts (RPKM, CPM, etc.).
//! *   Aggregation of signals from multiple BAM files.
//! *   Output to standard genomic formats.

use std::fmt::{Display, Formatter, Result as FmtResult};
use std::path::PathBuf;
use std::process::Command;

use ahash::HashMap;
use anyhow::{Context, Result};
use indicatif::ParallelProgressIterator;
use log::{error, info};
use noodles::{bam, core};
use polars::lazy::dsl::{col, cols, mean_horizontal, sum_horizontal};
use polars::prelude::*;
use rayon::prelude::*;
use rust_lapper::Lapper;
use tempfile;

use crate::bam_utils::{BamStats, Iv, get_bam_header, progress_bar};
use crate::genomic_intervals::IntervalMaker;
use crate::read_filter::BamReadFilter;
use crate::signal_normalization::NormalizationMethod;

/// Represents a single BAM file pileup settings.
pub struct BamPileup {
    file_path: PathBuf,
    bin_size: u64,
    norm_method: NormalizationMethod,
    scale_factor: f32,
    use_fragment: bool,
    filter: BamReadFilter,
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

                    // Merge adjacent bins if their count is equal.
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

/// Represents methods for aggregating signals across multiple BAM files.
pub enum AggregationMethod {
    None,
    Sum,
    Mean,
}

/// Represents a pileup from multiple BAM files.
///
/// The pileup is generated for each genomic region and then aggregated across the
/// provided BAM files using the chosen aggregation method.
pub struct MultiBamPileup {
    file_paths: Vec<PathBuf>,
    bin_size: u64,
    norm_method: NormalizationMethod,
    scale_factor: f32,
    use_fragment: bool,
    filters: Vec<BamReadFilter>,
    aggregation_method: AggregationMethod,
}

impl MultiBamPileup {
    /// Create a new [`MultiBamPileup`].
    pub fn new(
        file_paths: Vec<PathBuf>,
        bin_size: u64,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        use_fragment: bool,
        filters: Vec<BamReadFilter>,
        aggregation_method: AggregationMethod,
    ) -> Self {
        Self {
            file_paths,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
            filters,
            aggregation_method,
        }
    }

    /// Generate a DataFrame with a column for each BAM file's score per bin.
    pub fn pileup(&self) -> Result<DataFrame> {
        // Use the first file to obtain genome chunks and chromosome sizes.
        let bam_stats = BamStats::new(self.file_paths[0].clone())?;
        let genomic_chunks = bam_stats.genome_chunks(self.bin_size)?;
        let chromsizes_refid = bam_stats.chromsizes_ref_id()?;
        let n_total_chunks = genomic_chunks.len();

        info!("MultiBamPileup settings:");
        info!("BAM files: {:?}", self.file_paths);
        info!("Bin size: {}", self.bin_size);
        info!("Processing {} genomic chunks", n_total_chunks);
        info!("Normalization method: {:?}", self.norm_method);
        info!("Scaling factor: {}", self.scale_factor);
        info!("Using fragment for counting: {}", self.use_fragment);
        info!("Filtering using: {:?}", self.filters);

        // Process each genomic chunk in parallel.
        let pileup = genomic_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                n_total_chunks as u64,
                "Performing multi-BAM pileup".to_string(),
            ))
            .map(|region| {
                // For each region, process all BAM files in parallel.
                let pileup_vec: Vec<Vec<Iv>> = self
                    .file_paths
                    .par_iter()
                    .enumerate()
                    .map(|(i, bam_file)| {
                        pileup_chunk(
                            region.clone(),
                            bam_file.clone(),
                            self.bin_size,
                            chromsizes_refid.clone(),
                            self.use_fragment,
                            &self.filters[i],
                            false,
                        )
                    })
                    .collect();

                // Use the intervals from the first BAM file as the reference.
                let mut chrom = Vec::new();
                let mut start = Vec::new();
                let mut end = Vec::new();

                for iv in &pileup_vec[0] {
                    chrom.push(region.name().to_owned().to_string());
                    start.push(iv.start as u64);
                    end.push(iv.stop as u64);
                }

                // Build a DataFrame with one score column per BAM file.
                let mut df = DataFrame::new(vec![
                    Column::new("chrom".into(), chrom),
                    Column::new("start".into(), start),
                    Column::new("end".into(), end),
                ])?;

                for (i, ivs) in pileup_vec.into_iter().enumerate() {
                    let scores: Vec<u32> = ivs.into_iter().map(|iv| iv.val).collect();
                    let col_name = format!("score_{}", i);
                    df.with_column(Series::new(col_name.into(), scores))?;
                }

                Ok(df)
            })
            // Combine the DataFrames from all genomic chunks.
            .reduce(
                || Ok(DataFrame::empty()),
                |acc, result| match (acc, result) {
                    (Ok(mut acc_df), Ok(df)) => {
                        acc_df.vstack_mut(&df).context("Error stacking DataFrames")
                    }
                    (Err(e), _) | (_, Err(e)) => Err(e),
                },
            )?;

        Ok(pileup)
    }

    /// Normalize each BAM file's score column in the DataFrame.
    fn pileup_normalised(&self) -> Result<DataFrame> {
        let mut df = self.pileup()?;
        let score_cols: Vec<String> = (0..self.file_paths.len())
            .map(|i| format!("score_{}", i))
            .collect();

        for col_name in score_cols.iter() {
            let total: u64 = df
                .column(col_name)?
                .sum_reduce()?
                .as_any_value()
                .try_extract()?;
            let norm_factor =
                self.norm_method
                    .scale_factor(self.scale_factor, self.bin_size, total);
            let norm_scores = df.column(col_name)?.u32()? * norm_factor;
            df.replace(col_name, norm_scores)?;
        }
        Ok(df)
    }

    /// Aggregate the normalized scores from all BAM files using the selected method.
    fn pileup_aggregated(&self) -> Result<DataFrame> {
        let mut df = self.pileup_normalised()?;

        // Score column names are just numbers with score_i format.
        // Make a vector with these names to use in the aggregation.
        let score_cols: Vec<_> = (0..self.file_paths.len())
            .map(|i| format!("score_{}", i))
            // .map(|s| col(s))
            .collect();

        // Aggregate the score columns using the selected method.
        let df = match self.aggregation_method {
            AggregationMethod::Sum => df.with_column(
                df.select(&score_cols)?
                    .sum_horizontal(NullStrategy::Ignore)?
                    .context("Error summing scores")?
                    .with_name("score".into()),
            ),
            AggregationMethod::Mean => df.with_column(
                df.select(&score_cols)?
                    .mean_horizontal(NullStrategy::Ignore)?
                    .context("Error averaging scores")?
                    .with_name("score".into()),
            ),
        }?;
        // Remove the original per-file score columns.
        let df = df.drop_many(&score_cols);
        Ok(df)
    }

    /// Write the aggregated, normalized pileup as a bedGraph file.
    pub fn to_bedgraph(&self, outfile: PathBuf) -> Result<()> {
        info!("Writing bedGraph file to {}", outfile.display());
        let df = self.pileup_aggregated()?;
        let collapsed_df = collapse_equal_bins(df)?;
        write_bedgraph(collapsed_df, outfile)
    }

    /// Write the aggregated, normalized pileup as a BigWig file.
    pub fn to_bigwig(&self, outfile: PathBuf) -> Result<()> {
        let bam_stats = BamStats::new(self.file_paths[0].clone())?;
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
fn collapse_equal_bins(mut df: DataFrame) -> Result<DataFrame> {
    df = df.sort(["chrom", "start"], Default::default())?;
    let shifted_score = df.column("score")?.shift(1);
    let shifted_chrom = df.column("chrom")?.shift(1);
    let same_chrom = df.column("chrom")?.equal(&shifted_chrom)?;
    let same_score = df.column("score")?.equal(&shifted_score)?;
    let same_chrom_and_score = (same_chrom & same_score).cast(&DataType::Int8)?;

    // Compute a cumulative sum to define groups of identical rows.
    let group = cum_sum(&same_chrom_and_score, false)?
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

/// Generates a pileup for a given genomic region from a BAM file.
///
/// This function opens the specified BAM file, queries the given region,
/// creates intervals from the reads (using the provided filtering criteria),
/// and then aggregates the intervals into bins of size `bin_size` using a lapper.
/// If `collapse_equal_bins` is set to true, adjacent bins with the same count are merged.
///
/// # Arguments
///
/// * `region` - The genomic region for which the pileup is generated.
/// * `file_path` - The path to the BAM file.
/// * `bin_size` - The size of the bins to use for counting.
/// * `chromsizes_refid` - A mapping from reference IDs to chromosome sizes.
/// * `use_fragment` - Whether to count fragments instead of individual reads.
/// * `filter` - The read filter used to select which reads to count.
/// * `collapse_equal_bins` - If true, adjacent bins with identical counts are merged.
///
/// # Returns
///
/// A vector of `Iv` structs representing the pileup counts for each bin.
fn pileup_chunk(
    region: Region,
    file_path: PathBuf,
    bin_size: u64,
    chromsizes_refid: HashMap<usize, u64>,
    use_fragment: bool,
    filter: &BamReadFilter,
    collapse_equal_bins: bool,
) -> Result<Vec<Iv>> {
    // Open the BAM file and retrieve its header.
    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(file_path.clone())
        .context("Error opening BAM file")?;
    let header = get_bam_header(file_path.clone()).context("Error getting BAM header")?;

    // Query the BAM file for records overlapping the region.
    let records = reader
        .query(&header, &region)
        .context("Error querying BAM file")?;

    // Create intervals from valid records that pass filtering.
    let intervals: Vec<Iv> = records
        .filter_map(Result::ok)
        .filter_map(|record| {
            IntervalMaker::new(
                record,
                &header,
                &chromsizes_refid,
                filter,
                use_fragment,
                None,
                None,
            )
            .coords()
        })
        .map(|(start, end)| Iv {
            start,
            stop: end,
            val: 1,
        })
        .collect();

    // Build a lapper from the intervals.
    let lapper = Lapper::new(intervals);

    // Obtain the start and end coordinates for the region.
    let region_interval = region.interval();
    let region_start = region_interval
        .start()
        .context("Error getting interval start")?
        .get();
    let region_end = region_interval
        .end()
        .context("Error getting interval end")?
        .get();

    let mut bin_counts = Vec::new();
    let mut start = region_start;

    // Iterate over the region in bins.
    while start < region_end {
        let end = (start + bin_size).min(region_end);
        let count = lapper.count(start as usize, end as usize) as u32;

        if collapse_equal_bins {
            if let Some(last) = bin_counts.last_mut() {
                if last.val == count {
                    last.stop = end;
                } else {
                    bin_counts.push(Iv {
                        start,
                        stop: end,
                        val: count,
                    });
                }
            } else {
                bin_counts.push(Iv {
                    start,
                    stop: end,
                    val: count,
                });
            }
        } else {
            bin_counts.push(Iv {
                start,
                stop: end,
                val: count,
            });
        }
        start = end;
    }

    Ok(bin_counts)
}
