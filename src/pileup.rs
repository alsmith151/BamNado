//! This module implements pileup generation from one or multiple BAM files
//! and converts the resulting signal into bedGraph and BigWig formats. The
//! implementation is parallelized (using Rayon) and uses several libraries to
//! process genomic intervals and to normalize and aggregate counts.
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::ops::Deref;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use ahash::HashMap;
use anyhow::{Context, Result};
use indicatif::ParallelProgressIterator;
use log::{debug, error, info, warn};
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
    df = df.sort(["chrom", "start"], Default::default())
        .context("Failed to sort DataFrame for bedGraph output")?;
    
    let mut file = std::fs::File::create(&outfile)
        .context(format!("Failed to create output file: {}", outfile.display()))?;
    
    CsvWriter::new(&mut file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(&mut df)
        .context("Failed to write bedGraph file")?;
    
    Ok(())
}

/// Collapse adjacent bins that have equal chromosome and score values.
///
/// This function reduces the number of rows by merging adjacent bins with the
/// same score, which can make the output bedGraph file smaller.
fn collapse_equal_bins(
    df: DataFrame,
    score_columns: Option<Vec<PlSmallStr>>,
) -> Result<DataFrame> {
    let mut df = df.sort(["chrom", "start"], Default::default())
        .context("Failed to sort DataFrame for collapsing bins")?;
    
    let shifted_chromosome = df.column("chrom")?.shift(1);
    let same_chrom = df.column("chrom")?
        .equal(&shifted_chromosome)
        .context("Failed to compare chromosome columns")?;

    // For each row check that all the scores are the same in the shifted score columns as the original score columns
    let score_columns = score_columns.unwrap_or_else(|| vec!["score".into()]);

    let scores = df.select(score_columns.iter().map(|s| s.clone()))
        .context("Failed to select score columns")?;
    let scores_shifted = df.select(score_columns.iter().map(|s| &**s))
        .context("Failed to select score columns for shifting")?
        .shift(1);
    
    let mut same_scores = Vec::with_capacity(scores.width());
    for (score, score_shifted) in scores.iter().zip(scores_shifted.iter()) {
        let same_score = score.equal(score_shifted)
            .context("Failed to compare score columns")?;
        same_scores.push(same_score);
    }

    // Sum the boolean columns to get a single boolean column
    let same_chrom_and_score = same_scores
        .iter()
        .fold(same_chrom, |acc, x| acc & x.clone());

    // Compute a cumulative sum to define groups of identical rows.
    let group = cum_sum(&same_chrom_and_score.into_series(), false)
        .context("Failed to compute cumulative sum")?
        .with_name("groups".into())
        .into_column();

    let df = df
        .with_column(group)
        .context("Failed to add group column")?
        .clone()
        .lazy()
        .group_by(["groups"])
        .agg(&[
            col("chrom").first(),
            col("start").min(),
            col("end").max(),
            // Use appropriate aggregation based on score columns
            // If there are multiple score columns, keep their original values
            // as they should be the same within each group
            col("score").first(),
        ])
        .collect()
        .context("Failed to group and aggregate data")?;

    Ok(df)
}

/// Convert a bedGraph file to BigWig format using the external command
/// `bedGraphToBigWig`.
fn convert_bedgraph_to_bigwig(
    bedgraph_path: &std::path::Path,
    chromsizes_path: &std::path::Path,
    outfile: &PathBuf,
) -> Result<()> {
    // Use Command builder pattern for better readability
    let output = Command::new("bedGraphToBigWig")
        .arg(bedgraph_path)
        .arg(chromsizes_path)
        .arg(outfile)
        .output()
        .context(format!(
            "Failed to execute bedGraphToBigWig on file: {}",
            bedgraph_path.display()
        ))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        error!("Error converting bedGraph to BigWig:");
        error!("{}", stderr);
        anyhow::bail!(
            "Conversion to BigWig failed: {}",
            stderr
        );
    } else {
        info!("BigWig file successfully written to {}", outfile.display());
    }
    Ok(())
}

/// Represents a configuration for generating pileup data from a BAM file.
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
    // Counter for failed chunks during parallel processing
    failed_chunks: AtomicUsize,
}

impl Display for BamPileup {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        writeln!(f, "Pileup Settings:")?;
        writeln!(f, "BAM file: {}", self.file_path.display())?;
        writeln!(f, "Bin size: {}", self.bin_size)?;
        writeln!(f, "Normalization method: {:?}", self.norm_method)?;
        writeln!(f, "Scaling factor: {}", self.scale_factor)?;
        writeln!(f, "Using fragment for counting: {}", self.use_fragment)?;
        writeln!(f, "Collapse equal bins: {}", self.collapse)?;
        write!(f, "Filtering using: \n{}\n", self.filter)
    }
}

impl BamPileup {
    /// Create a new [`BamPileup`].
    ///
    /// # Arguments
    ///
    /// * `file_path` - Path to the BAM file
    /// * `bin_size` - Size of bins for counting reads
    /// * `norm_method` - Method to normalize counts
    /// * `scale_factor` - Factor to scale normalized counts
    /// * `use_fragment` - Whether to use fragments instead of reads
    /// * `filter` - Filter for reads to include
    /// * `collapse` - Whether to collapse adjacent bins with equal scores
    ///
    /// # Returns
    ///
    /// A new [`BamPileup`] instance if validation passes
    pub fn new(
        file_path: PathBuf,
        bin_size: u64,
        norm_method: NormalizationMethod,
        scale_factor: f32,
        use_fragment: bool,
        filter: BamReadFilter,
        collapse: bool, // Fixed typo from 'colapse' to 'collapse'
    ) -> Result<Self> {
        // Validate inputs
        if bin_size == 0 {
            return Err(anyhow::anyhow!("Bin size must be greater than 0"));
        }
        
        if scale_factor <= 0.0 {
            return Err(anyhow::anyhow!("Scale factor must be greater than 0"));
        }
        
        // Validate that the BAM file exists and is readable
        if !file_path.exists() {
            return Err(anyhow::anyhow!("BAM file does not exist: {}", file_path.display()));
        }
        
        // Check if the corresponding .bai index file exists
        let index_path = PathBuf::from(format!("{}.bai", file_path.display()));
        if !index_path.exists() {
            warn!("BAM index file not found: {}. This may cause issues with random access.", index_path.display());
        }
        
        Ok(Self {
            file_path,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
            filter,
            collapse,
            failed_chunks: AtomicUsize::new(0),
        })
    }

    /// Generate pileup intervals for the BAM file.
    ///
    /// Returns a map from chromosome names to a vector of intervals (with
    /// start, stop and count) for each genomic chunk.
    fn pileup(&self) -> Result<HashMap<String, Vec<Iv>>> {
        let bam_stats = BamStats::new(self.file_path.clone())
            .context(format!("Failed to get BAM stats for file: {}", self.file_path.display()))?;
            
        let genomic_chunks = bam_stats.genome_chunks(self.bin_size)
            .context("Failed to generate genomic chunks")?;
            
        let chromsizes_refid = bam_stats.chromsizes_ref_id()
            .context("Failed to get chromosome sizes with reference IDs")?;
            
        let n_total_chunks = genomic_chunks.len();

        info!("{}", self);
        info!("Processing {} genomic chunks", n_total_chunks);

        let header = get_bam_header(&self.file_path)
            .context(format!("Failed to get BAM header from file: {}", self.file_path.display()))?;

        // Reset failed chunks counter
        self.failed_chunks.store(0, Ordering::SeqCst);

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
                    .build_from_path(&self.file_path)
                    .context(format!(
                        "Failed to open BAM file: {}",
                        self.file_path.display()
                    ))?;

                // Query for reads overlapping the region.
                let records = reader.query(&header, &region)
                    .context(format!("Failed to query region: {:?}", region))?;

                // Create intervals from each read that passes filtering.
                let intervals: Vec<Iv> = records
                    .into_iter()
                    .filter_map(|record| {
                        match record {
                            Ok(r) => Some(r),
                            Err(e) => {
                                warn!("Failed to read record: {}", e);
                                None
                            }
                        }
                    })
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

                    let new_interval = rust_lapper::Interval {
                        start,
                        stop: end,
                        val: count as u32,
                    };
                    
                    if self.collapse && !bin_counts.is_empty() {
                        let last = bin_counts.last_mut().unwrap();
                        if last.val == count as u32 {
                            last.stop = end; // Extend the last interval
                        } else {
                            bin_counts.push(new_interval);
                        }
                    } else {
                        bin_counts.push(new_interval);
                    }
                    
                    start = end;
                }

                Ok((region.name().to_owned().to_string(), bin_counts))
            })
            // Combine the results from parallel threads.
            .fold(
                || HashMap::default(),
                |mut acc, result: Result<(String, Vec<Iv>)>| {
                    match result {
                        Ok((chrom, intervals)) => {
                            acc.entry(chrom).or_insert_with(Vec::new).extend(intervals);
                        }
                        Err(e) => {
                            error!("Error processing chunk: {}", e);
                            self.failed_chunks.fetch_add(1, Ordering::SeqCst);
                        }
                    }
                    acc
                },
            )
            .reduce(
                || HashMap::default(),
                |mut acc, map| {
                    for (key, value) in map {
                        acc.entry(key)
                          .or_insert_with(Vec::new)
                          .extend(value);
                    }
                    acc
                },
            );

        let failed_count = self.failed_chunks.load(Ordering::SeqCst);
        if failed_count > 0 {
            warn!("{} genomic chunks failed to process", failed_count);
        }

        info!("Pileup complete");
        info!("Read filtering statistics: {}", self.filter.stats());

        Ok(pileup)
    }

    /// Convert the pileup intervals into a Polars DataFrame.
    ///
    /// The DataFrame will have columns "chrom", "start", "end" and "score".
    fn pileup_to_polars(&self) -> Result<DataFrame> {
        let pileup = self.pileup()
            .context("Failed to generate pileup")?;

        // Calculate total capacity needed
        let total_size: usize = pileup.values().map(|v| v.len()).sum();
        
        // Pre-allocate vectors with exact capacity
        let mut chroms = Vec::with_capacity(total_size);
        let mut starts = Vec::with_capacity(total_size);
        let mut ends = Vec::with_capacity(total_size);
        let mut scores = Vec::with_capacity(total_size);
        
        // Process each chromosome and flatten into the vectors
        for (chrom, intervals) in pileup {
            for iv in intervals {
                chroms.push(chrom.clone());
                starts.push(iv.start as u64);
                ends.push(iv.stop as u64);
                scores.push(iv.val as u32);
            }
        }
        
        // Create the DataFrame using Series for better performance
        let df = DataFrame::new(vec![
            Series::new("chrom".into(), chroms).into(),
            Series::new("start".into(), starts).into(),
            Series::new("end".into(), ends).into(),
            Series::new("score".into(), scores).into(),
        ]).context("Failed to create pileup DataFrame")?;
        
        Ok(df)
    }

    /// Normalize the pileup signal using the provided normalization method.
    fn pileup_normalised(&self) -> Result<DataFrame> {
        let mut df = self.pileup_to_polars()
            .context("Failed to convert pileup to DataFrame")?;
            
        // Get the total counts across all bins.
        let n_total_counts: u64 = df
            .column("score")?
            .sum_reduce()?
            .as_any_value()
            .try_extract()
            .context("Failed to extract total count")?;
            
        info!("Total counts: {}", n_total_counts);

        // Compute the normalization factor.
        let norm_factor =
            self.norm_method
                .scale_factor(self.scale_factor, self.bin_size, n_total_counts);
                
        // Multiply the score column by the normalization factor.
        let norm_scores = df.column("score")?
            .u32()
            .context("Failed to convert score column to u32")?
            * norm_factor;
            
        df.replace("score", norm_scores)
            .context("Failed to replace score column with normalized values")?;
            
        Ok(df)
    }

    /// Write the normalized pileup as a bedGraph file.
    pub fn to_bedgraph(&self, outfile: PathBuf) -> Result<()> {
        info!("Writing bedGraph file to {}", outfile.display());
        let df = self.pileup_normalised()
            .context("Failed to normalize pileup")?;
            
        write_bedgraph(df, outfile)
            .context("Failed to write bedGraph file")?;
            
        Ok(())
    }

    /// Write the normalized pileup as a BigWig file.
    ///
    /// This function writes a temporary bedGraph file and converts it to BigWig.
    pub fn to_bigwig(&self, outfile: PathBuf) -> Result<()> {
        let bam_stats = BamStats::new(self.file_path.clone())
            .context(format!("Failed to get BAM stats for file: {}", self.file_path.display()))?;
            
        let chromsizes_file = tempfile::NamedTempFile::new()
            .context("Failed to create temporary file for chromosome sizes")?;
            
        let chromsizes_path = chromsizes_file.path();
        bam_stats.write_chromsizes(chromsizes_path.to_path_buf())
            .context("Failed to write chromosome sizes to temporary file")?;

        let bedgraph_file = tempfile::NamedTempFile::new()
            .context("Failed to create temporary file for bedGraph")?;
            
        let bedgraph_path = bedgraph_file.path();
        self.to_bedgraph(bedgraph_path.to_path_buf())
            .context("Failed to write to temporary bedGraph file")?;

        info!("Converting bedGraph to BigWig file");
        convert_bedgraph_to_bigwig(bedgraph_path, chromsizes_path, &outfile)
            .context("Failed to convert bedGraph to BigWig")?;
            
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
    /// * `pileups` - Vector of BamPileup instances to process in parallel
    ///
    /// # Returns
    ///
    /// A new [`MultiBamPileup`] instance if validation passes
    pub fn new(pileups: Vec<BamPileup>) -> Result<Self> {
        let instance = Self {
            bam_pileups: pileups,
        };
        
        // Validate on construction for early failure
        instance.validate_bam_pileups()?;
        
        Ok(instance)
    }

    /// Validate that all BAM pileups have compatible settings.
    fn validate_bam_pileups(&self) -> Result<()> {
        if self.bam_pileups.is_empty() {
            return Err(anyhow::anyhow!("No BAM pileups provided"));
        }
        
        // Check that all BAM pileups have the same bin size.
        let bin_sizes: Vec<u64> = self.bam_pileups.iter().map(|p| p.bin_size).collect();
        if !bin_sizes.iter().all(|&x| x == bin_sizes[0]) {
            return Err(anyhow::anyhow!(
                "All BAM pileups must have the same bin size, found: {:?}", 
                bin_sizes
            ));
        } else {
            debug!("All BAM pileups have the same bin size: {}", bin_sizes[0]);
        }

        // Check that all BAM pileups have the collapse option set to false
        let collapse: Vec<bool> = self.bam_pileups.iter().map(|p| p.collapse).collect();
        if !collapse.iter().all(|x| *x == false) {
            return Err(anyhow::anyhow!(
                "All BAM pileups must have collapse set to false for multi-BAM processing"
            ));
        } else {
            debug!("All BAM pileups have collapse set to false");
        }

        Ok(())
    }

    /// Generate pileup intervals for all BAM files.
    ///
    /// Returns a DataFrame with columns for chromosome, start, end,
    /// and scores for each BAM file.
    fn pileup(&self) -> Result<DataFrame> {
        let mut dfs = Vec::with_capacity(self.bam_pileups.len());
        let names: Vec<PlSmallStr> = self
            .bam_pileups
            .iter()
            .map(|p| {
                p.file_path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or_else(|| p.file_path.to_str().unwrap_or("unknown"))
                    .into()
            })
            .collect::<Vec<_>>();

        // Generate pileups for each BAM file
        for (i, pileup) in self.bam_pileups.iter().enumerate() {
            let df = pileup.pileup_to_polars()
                .context(format!("Failed to generate pileup for BAM file: {}", pileup.file_path.display()))?;
            dfs.push(df);
        }

        // Initialize with the first DataFrame and rename its score column
        let df = dfs[0]
            .clone()
            .lazy()
            .rename(["score"], [names[0].clone()], true);
            
        // Extract score columns from remaining DataFrames
        let scores = dfs
            .iter()
            .skip(1)
            .zip(names.iter().skip(1))
            .map(|(df, n)| {
                df.column("score")
                  .map(|c| c.clone().with_name(n.clone()))
                  .context(format!("Failed to get score column for {}", n))
            })
            .collect::<Result<Vec<Column>>>()?;
            
        // Create a DataFrame with just the score columns
        let scores_df = DataFrame::new(scores)
            .context("Failed to create DataFrame from score columns")?
            .lazy();
            
        // Join all DataFrames
        let df = df.join(
            scores_df,
            [col("chrom"), col("start"), col("end")],
            [col("chrom"), col("start"), col("end")],
            JoinArgs::new(JoinType::Full),
        )
        .collect()
        .context("Failed to join DataFrames")?;
        
        // Collapse equal bins if requested
        let df = collapse_equal_bins(df, Some(names))
            .context("Failed to collapse equal bins")?;

        Ok(df)
    }

    /// Write the multi-BAM pileup to a TSV file.
    pub fn to_tsv(&self, outfile: &PathBuf) -> Result<()> {
        let mut df = self.pileup()
            .context("Failed to generate multi-BAM pileup")?;
            
        // Normalize scores if needed
        for pileup in &self.bam_pileups {
            // Implementation for normalizing individual columns would go here
            // This is currently missing in the original code
        }
            
        let mut file = std::fs::File::create(outfile)
            .context(format!("Failed to create output file: {}", outfile.display()))?;
            
        CsvWriter::new(&mut file)
            .include_header(true)
            .with_separator(b'\t')
            .finish(&mut df)
            .context("Failed to write TSV file")?;

        info!("Multi-BAM pileup written to {}", outfile.display());
        Ok(())
    }
    
    /// Write the multi-BAM pileup as a BigWig file.
    /// 
    /// This would be a new method to support BigWig output for multi-BAM pileups.
    pub fn to_bigwig(&self, outfile: &PathBuf) -> Result<()> {
        // Implementation would be similar to BamPileup::to_bigwig
        // but handling multiple score columns
        unimplemented!("BigWig output for multi-BAM pileups is not yet implemented");
    }
}