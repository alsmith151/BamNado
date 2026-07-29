//! # Bin Counts Module
//!
//! This module computes a shared genomic bin grid across one or more BAM files and
//! counts filtered reads (or fragments) per bin per sample. It supports:
//! *   A single shared bin grid across samples, required for between-sample comparison.
//! *   Parallel per-chunk counting, reusing the same binning machinery as `coverage_analysis`.
//! *   Conversion of the resulting count matrix to a Polars `DataFrame` for QC / M-A plots.

use std::path::PathBuf;

use crate::bam_utils::{BamStats, Iv, bin_intervals, get_bam_header, is_scaffold, progress_bar};
use crate::genomic_intervals::IntervalMaker;
use crate::read_filter::BamReadFilter;
use ahash::HashMap;
use anyhow::{Context, Result, bail};
use indicatif::ParallelProgressIterator;
use ndarray::Array2;
use noodles::bam;
use polars::prelude::*;
use rayon::prelude::*;
use rust_lapper::Lapper;

/// One genomic bin in the count matrix.
///
/// Coordinates are 1-based (noodles/SAM convention), matching `BamStats`.
#[derive(Debug, Clone)]
pub struct BinRegion {
    /// Chromosome / reference sequence name.
    pub chrom: String,
    /// 1-based inclusive start position.
    pub start: u64,
    /// 1-based inclusive end position.
    pub end: u64,
}

/// A sample × bin read-count matrix, produced by [`count_bins`].
pub struct BinCounts {
    /// Sample names, in column order matching `counts`.
    pub sample_names: Vec<String>,
    /// The bin size used to build the grid.
    pub bin_size: u64,
    /// The bin grid, in row order matching `counts`.
    pub bins: Vec<BinRegion>,
    /// Read counts, shape `(n_bins, n_samples)`.
    pub counts: Array2<u64>,
    /// Post-filter library size (total retained reads/fragments) per sample.
    pub library_sizes: Vec<u64>,
}

impl BinCounts {
    /// Convert the bin count matrix into a Polars `DataFrame` for QC (e.g. M-A plots).
    ///
    /// # Returns
    ///
    /// A `DataFrame` with columns `chrom`, `start`, `end`, and one `f64` column per sample,
    /// named after the corresponding entry in `sample_names`.
    pub fn to_polars(&self) -> Result<DataFrame> {
        let n_bins = self.bins.len();
        let mut chroms = Vec::with_capacity(n_bins);
        let mut starts = Vec::with_capacity(n_bins);
        let mut ends = Vec::with_capacity(n_bins);
        for bin in &self.bins {
            chroms.push(bin.chrom.clone());
            starts.push(bin.start);
            ends.push(bin.end);
        }

        let mut columns = vec![
            Column::new("chrom".into(), chroms),
            Column::new("start".into(), starts),
            Column::new("end".into(), ends),
        ];

        for (s, name) in self.sample_names.iter().enumerate() {
            let col: Vec<f64> = self.counts.column(s).iter().map(|&c| c as f64).collect();
            columns.push(Column::new(name.as_str().into(), col));
        }

        Ok(DataFrame::new(columns)?)
    }
}

/// Bin-count reads (or fragments) across one or more BAM files onto a shared genomic bin grid.
///
/// All input BAM files must share an identical chromosome name → length mapping; this is
/// required so that a single bin grid can be used to compare samples bin-for-bin.
///
/// # Arguments
///
/// * `bam_paths` - Input BAM files, one per sample.
/// * `sample_names` - Sample names, same length and order as `bam_paths`.
/// * `bin_size` - Genomic bin width in base pairs.
/// * `filters` - Per-sample read filter, same length and order as `bam_paths`.
/// * `use_fragment` - Count fragments (read pairs) instead of individual read alignments.
/// * `ignore_scaffold_chromosomes` - Skip scaffold / unplaced chromosomes when building the grid.
///
/// # Returns
///
/// A [`BinCounts`] with the count matrix and per-sample post-filter library sizes.
pub fn count_bins(
    bam_paths: &[PathBuf],
    sample_names: &[String],
    bin_size: u64,
    filters: &[BamReadFilter],
    use_fragment: bool,
    ignore_scaffold_chromosomes: bool,
) -> Result<BinCounts> {
    anyhow::ensure!(
        !bam_paths.is_empty(),
        "At least one BAM file is required for bin counting"
    );
    anyhow::ensure!(
        bam_paths.len() == sample_names.len(),
        "Number of BAM paths ({}) does not match number of sample names ({})",
        bam_paths.len(),
        sample_names.len()
    );
    anyhow::ensure!(
        bam_paths.len() == filters.len(),
        "Number of BAM paths ({}) does not match number of filters ({})",
        bam_paths.len(),
        filters.len()
    );

    // Build per-sample BamStats and validate a shared chromosome grid.
    let bam_stats: Vec<BamStats> = bam_paths
        .iter()
        .map(|p| {
            BamStats::new(p.clone())
                .with_context(|| format!("Failed to read BAM stats for {}", p.display()))
        })
        .collect::<Result<Vec<_>>>()?;

    let reference_chromsizes = bam_stats[0].chromsizes_ref_name()?;
    for (stats, path) in bam_stats.iter().zip(bam_paths.iter()).skip(1) {
        let chromsizes = stats.chromsizes_ref_name()?;
        if chromsizes != reference_chromsizes {
            bail!(
                "Chromosome name/length mismatch between BAM files: {} does not match {}. \
                 bam-normalize requires all input BAM files to share an identical set of \
                 reference sequences so that a common bin grid can be used.",
                path.display(),
                bam_paths[0].display()
            );
        }
    }

    // Build the global bin grid: one contiguous block of bins per chromosome, in
    // deterministic (sorted) chromosome name order.
    let mut chrom_names: Vec<String> = reference_chromsizes.keys().cloned().collect();
    if ignore_scaffold_chromosomes {
        chrom_names.retain(|name| !is_scaffold(name));
    }
    chrom_names.sort();

    let mut bins = Vec::new();
    let mut chrom_offsets: HashMap<String, usize> = HashMap::default();
    for chrom in &chrom_names {
        chrom_offsets.insert(chrom.clone(), bins.len());
        let length = reference_chromsizes[chrom];
        let mut start = 1u64;
        while start <= length {
            let end = (start + bin_size - 1).min(length);
            bins.push(BinRegion {
                chrom: chrom.clone(),
                start,
                end,
            });
            start = end + 1;
        }
    }
    let n_bins = bins.len();
    anyhow::ensure!(
        n_bins > 0,
        "No bins produced; check --bin-size and --ignore-scaffolds"
    );

    let mut counts = Array2::<u64>::zeros((n_bins, bam_paths.len()));
    let mut library_sizes = Vec::with_capacity(bam_paths.len());

    for (s, ((bam_path, stats), filter)) in bam_paths
        .iter()
        .zip(bam_stats.iter())
        .zip(filters.iter())
        .enumerate()
    {
        let header = get_bam_header(bam_path.clone())?;
        let chromsizes_refid = stats.chromsizes_ref_id()?;
        let genomic_chunks: Vec<_> = stats
            .genome_chunks(bin_size)?
            .into_iter()
            .filter(|region| {
                if !ignore_scaffold_chromosomes {
                    return true;
                }
                let name = region.name().to_owned().to_string();
                !is_scaffold(&name)
            })
            .collect();
        let n_total_chunks = genomic_chunks.len();

        let chunk_results: Vec<(String, usize, u32)> = genomic_chunks
            .into_par_iter()
            .progress_with(progress_bar(
                n_total_chunks as u64,
                format!("Counting bins for {}", sample_names[s]),
            ))
            .map(|region| -> Result<Vec<(String, usize, u32)>> {
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(bam_path)
                    .context("Failed to open BAM file")?;
                let records = reader.query(&header, &region)?;
                let intervals: Vec<Iv> = records
                    .records()
                    .filter_map(|record| record.ok())
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
                        .map(|(start_pos, end_pos, _dtlen)| Iv {
                            start: start_pos,
                            stop: end_pos,
                            val: 1,
                        })
                    })
                    .collect();

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

                let chrom_name = region.name().to_owned().to_string();
                let bin_counts = bin_intervals(&lapper, region_start, region_end, bin_size, false);
                Ok(bin_counts
                    .into_iter()
                    .map(|iv| (chrom_name.clone(), iv.start, iv.val))
                    .collect())
            })
            .fold(
                || Ok(Vec::new()),
                |acc: Result<Vec<(String, usize, u32)>>,
                 result: Result<Vec<(String, usize, u32)>>| {
                    let mut acc = acc?;
                    acc.extend(result?);
                    Ok(acc)
                },
            )
            .reduce(
                || Ok(Vec::new()),
                |acc: Result<Vec<(String, usize, u32)>>, b: Result<Vec<(String, usize, u32)>>| {
                    let mut acc = acc?;
                    acc.extend(b?);
                    Ok(acc)
                },
            )?;

        for (chrom, start, count) in chunk_results {
            let offset = *chrom_offsets.get(&chrom).ok_or_else(|| {
                anyhow::anyhow!("Chromosome {chrom} not part of the shared bin grid")
            })?;
            let bin_idx_in_chrom = ((start as u64 - 1) / bin_size) as usize;
            counts[[offset + bin_idx_in_chrom, s]] = count as u64;
        }

        library_sizes.push(filter.stats().n_reads_after_filtering());
    }

    Ok(BinCounts {
        sample_names: sample_names.to_vec(),
        bin_size,
        bins,
        counts,
        library_sizes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_filter::BamReadFilter;

    fn get_test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("Failed to get workspace root")
            .join("test/data")
    }

    #[test]
    fn test_count_bins_same_bam_twice_gives_identical_columns() {
        let bam_file = get_test_data_dir().join("test.bam");
        let bam_paths = vec![bam_file.clone(), bam_file];
        let sample_names = vec!["a".to_string(), "b".to_string()];
        let filters = vec![BamReadFilter::default(), BamReadFilter::default()];

        let counts = count_bins(&bam_paths, &sample_names, 10_000, &filters, false, true)
            .expect("count_bins failed");

        assert!(!counts.bins.is_empty(), "Expected at least one bin");
        assert_eq!(counts.library_sizes[0], counts.library_sizes[1]);
        assert!(counts.library_sizes[0] > 0);

        let total: u64 = counts.counts.iter().sum();
        assert!(total > 0, "Expected non-zero counts");

        for row in counts.counts.rows() {
            assert_eq!(row[0], row[1]);
        }
    }

    #[test]
    fn test_count_bins_mismatched_bam_paths_and_sample_names_errors() {
        let bam_file = get_test_data_dir().join("test.bam");
        let bam_paths = vec![bam_file];
        let sample_names = vec!["a".to_string(), "b".to_string()];
        let filters = vec![BamReadFilter::default()];

        let result = count_bins(&bam_paths, &sample_names, 10_000, &filters, false, true);
        assert!(result.is_err());
    }
}
