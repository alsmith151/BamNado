//!  BigWig file comparison utilities.
//! Aimed at comparing coverage BigWig files, e.g., from different samples or conditions.

use crate::bam_utils::progress_bar;
use anyhow::Result;
use bigtools::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};
use indicatif::ProgressIterator;
use itertools::Itertools;
use log::info;
use ndarray::prelude::*;
use polars::prelude::*;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;
use tempfile;

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum Comparison {
    Subtraction,
    Ratio,
    LogRatio,
}

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum AggregationMode {
    Sum,
    Mean,
    Median,
    Max,
    Min,
}

/// Read a range from a BigWig file into a provided array chunk.
/// # Arguments
/// * `bw` - Path to the BigWig file.
/// * `chrom` - Chromosome name.
/// * `start` - Start position (0-based).
/// * `end` - End position (0-based, exclusive).
/// * `array_chunk` - Mutable slice to fill with values.
/// # Returns
/// A Result indicating success or failure.
fn read_bp_range_into_array_chunk(
    bw: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    array_chunk: &mut [f32],
) -> Result<()> {
    let mut bw = bigtools::BigWigRead::open_file(bw)?;
    let values = bw.get_interval(chrom, start, end)?;
    for interval_res in values {
        let interval = interval_res?;
        // clip interval to requested region
        let s = std::cmp::max(interval.start, start) as usize;
        let e = std::cmp::min(interval.end, end) as usize;
        if s >= e {
            continue;
        }
        let base_offset = start as usize;
        // iterate only inside clipped interval
        for pos in s..e {
            array_chunk[pos - base_offset] = interval.value;
        }
    }
    Ok(())
}

/// Compare two BigWig files and report differences in as a new BigWig file.
/// # Arguments
/// * `bw1_path` - Path to the first BigWig file.
/// * `bw2_path` - Path to the second BigWig file.
/// * `output_path` - Path to the output BigWig file containing differences.
/// * `chunksize` - Optional chunk size for processing (in base pairs).
/// # Returns
/// A Result indicating success or failure.
pub fn compare_bigwigs<P>(
    bw1_path: &P,
    bw2_path: &Path,
    output_path: &Path,
    comparison: Comparison,
    bin_size: u32,
    chunksize: Option<usize>,
    pseudocount: Option<f64>,
) -> Result<()>
where
    P: AsRef<Path> + std::fmt::Debug + Send + Sync,
{
    // Initially need to figure out the best genomic chunk size to extract from both files.
    // On HPCs (most use cases) IO is the bottleneck, so we want to minimize the number of reads.
    // Initially get the chromosome sizes
    let chrom_info = bigtools::BigWigRead::open_file(bw1_path.as_ref())?
        .chroms()
        .to_owned()
        .into_iter()
        .sorted_by(|a, b| a.name.cmp(&b.name))
        .collect::<Vec<_>>();

    let chunk_size = chunksize.unwrap_or(1_000_000); // Default to 1Mb chunks if not specified

    let dfs: Vec<DataFrame> = chrom_info
        .iter()
        .progress_with(progress_bar(
            chrom_info.len() as u64,
            "Processing chromosomes".to_string(),
        ))
        .map(|chromosome| {
            let mut arr1 = Array1::<f32>::zeros(chromosome.length as usize);
            let mut arr2 = Array1::<f32>::zeros(chromosome.length as usize);

            arr1.as_slice_mut()
                .expect("contiguous")
                .par_chunks_mut(chunk_size)
                .zip(
                    arr2.as_slice_mut()
                        .expect("contiguous")
                        .par_chunks_mut(chunk_size),
                )
                .enumerate()
                .try_for_each(|(chunk_idx, (a1_chunk, a2_chunk))| -> Result<()> {
                    let start = chunk_idx * chunk_size;
                    let end = std::cmp::min(start + chunk_size, chromosome.length as usize);
                    // Fill a1_chunk and a2_chunk directly (bp resolution)
                    read_bp_range_into_array_chunk(
                        bw1_path.as_ref(),
                        &chromosome.name,
                        start as u32,
                        end as u32,
                        a1_chunk,
                    )?;
                    read_bp_range_into_array_chunk(
                        bw2_path,
                        &chromosome.name,
                        start as u32,
                        end as u32,
                        a2_chunk,
                    )?;

                    Ok(())
                })?;

            // Perform comparison
            let pc = pseudocount.unwrap_or(1e-12) as f32;
            let diff = match comparison {
                Comparison::Subtraction => &arr1 - &arr2,
                Comparison::Ratio => &arr1 / (&arr2 + pc),
                Comparison::LogRatio => ((&arr1 + pc) / (&arr2 + pc)).mapv(|x| x.ln()),
            };

            // Create a bedgraph style output for this chromosome as a polars DataFrame
            // Binning
            let n_bins = (chromosome.length as f32 / bin_size as f32).ceil() as usize;
            let mut chroms = Vec::with_capacity(n_bins);
            let mut starts = Vec::with_capacity(n_bins);
            let mut ends = Vec::with_capacity(n_bins);
            let mut values: Vec<f32> = Vec::with_capacity(n_bins);
            for bin_idx in 0..n_bins {
                let start = bin_idx * bin_size as usize;
                let end = std::cmp::min(start + bin_size as usize, chromosome.length as usize);
                chroms.push(chromosome.name.clone());
                starts.push(start as u32);
                ends.push(end as u32);
                values.push(diff.slice(s![start..end]).mean().unwrap_or(0.0)); // Handle empty slices
            }

            let mut df = df![
                "chrom" => chroms,
                "start" => starts,
                "end" => ends,
                "score" => values,
            ]?;
            df.sort_in_place(["chrom", "start"], SortMultipleOptions::default())?;
            Ok(df)
        })
        .collect::<Result<Vec<_>>>()?;

    // Concatenate all DataFrames
    if dfs.is_empty() {
        return Ok(());
    }
    let mut final_df = dfs[0].clone();
    for df in dfs.iter().skip(1) {
        final_df.vstack_mut(df)?;
    }
    // Write chrom sizes to temp file
    let chromsizes_file = tempfile::NamedTempFile::new()?;
    {
        let mut writer = std::io::BufWriter::new(&chromsizes_file);
        for chrom in &chrom_info {
            writeln!(writer, "{}\t{}", chrom.name, chrom.length)?;
        }
    }

    // Write bedgraph to temp file
    info!("Writing temporary bedgraph file");
    let bedgraph_file = tempfile::NamedTempFile::new()?;
    let bedgraph_path = bedgraph_file.path();

    // Use CsvWriter to write tab-separated values
    CsvWriter::new(&bedgraph_file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(&mut final_df)?;

    // Convert to BigWig
    let args = bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs {
        bedgraph: bedgraph_path.to_string_lossy().to_string(),
        chromsizes: chromsizes_file.path().to_string_lossy().to_string(),
        output: output_path.to_string_lossy().to_string(),
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
    info!("Writing output BigWig to {:?}", output_path);
    bigtools::utils::cli::bedgraphtobigwig::bedgraphtobigwig(args)
        .map_err(|e| anyhow::anyhow!("Error converting bedgraph to bigwig: {}", e))?;

    Ok(())
}

/// Aggregate multiple BigWig files into a single output BigWig.
/// # Arguments
/// * `bw_paths` - Vector of paths to BigWig files to aggregate.
/// * `output_path` - Path to the output BigWig file.
/// * `aggregation_mode` - The aggregation method (Sum, Mean, Median, Max, Min).
/// * `bin_size` - Size of bins for binning the output signal.
/// * `chunksize` - Optional chunk size for processing (in base pairs).
/// * `pseudocount` - Optional pseudocount value for aggregation (added to all values).
/// # Returns
/// A Result indicating success or failure.
pub fn aggregate_bigwigs<P>(
    bw_paths: &[P],
    output_path: &Path,
    aggregation_mode: AggregationMode,
    bin_size: u32,
    chunksize: Option<usize>,
    pseudocount: Option<f64>,
) -> Result<()>
where
    P: AsRef<Path> + std::fmt::Debug + Send + Sync,
{
    if bw_paths.is_empty() {
        return Err(anyhow::anyhow!("At least one BigWig file is required"));
    }

    // Get chromosome information from the first file
    let chrom_info = bigtools::BigWigRead::open_file(bw_paths[0].as_ref())?
        .chroms()
        .to_owned()
        .into_iter()
        .sorted_by(|a, b| a.name.cmp(&b.name))
        .collect::<Vec<_>>();

    let chunk_size = chunksize.unwrap_or(1_000_000); // Default to 1Mb chunks if not specified

    let dfs: Vec<DataFrame> = chrom_info
        .iter()
        .progress_with(progress_bar(
            chrom_info.len() as u64,
            "Processing chromosomes".to_string(),
        ))
        .map(|chromosome| {
            // Allocate arrays for all input BigWigs at bp resolution
            let mut arrays: Vec<Array1<f32>> = (0..bw_paths.len())
                .map(|_| Array1::<f32>::zeros(chromosome.length as usize))
                .collect();

            // Fill arrays in parallel chunks
            arrays
                .iter_mut()
                .zip(bw_paths.iter())
                .enumerate()
                .try_for_each(|(_bw_idx, (arr, bw_path))| {
                    arr.as_slice_mut()
                        .expect("contiguous")
                        .par_chunks_mut(chunk_size)
                        .enumerate()
                        .try_for_each(|(chunk_idx, chunk)| -> Result<()> {
                            let start = chunk_idx * chunk_size;
                            let end = std::cmp::min(start + chunk_size, chromosome.length as usize);
                            read_bp_range_into_array_chunk(
                                bw_path.as_ref(),
                                &chromosome.name,
                                start as u32,
                                end as u32,
                                chunk,
                            )?;
                            Ok(())
                        })?;
                    Ok::<(), anyhow::Error>(())
                })?;

            // Perform aggregation
            let pc = pseudocount.unwrap_or(0.0) as f32;

            // Add pseudocount to all arrays
            let arrays: Vec<Array1<f32>> = arrays.into_iter().map(|arr| &arr + pc).collect();

            let aggregated = match aggregation_mode {
                AggregationMode::Sum => arrays
                    .iter()
                    .skip(1)
                    .fold(arrays[0].clone(), |acc, arr| &acc + arr),
                AggregationMode::Mean => {
                    let sum = arrays
                        .iter()
                        .skip(1)
                        .fold(arrays[0].clone(), |acc, arr| &acc + arr);
                    sum / arrays.len() as f32
                }
                AggregationMode::Max => {
                    arrays.iter().skip(1).fold(arrays[0].clone(), |acc, arr| {
                        acc.iter().zip(arr.iter()).map(|(a, b)| a.max(*b)).collect()
                    })
                }
                AggregationMode::Min => {
                    arrays.iter().skip(1).fold(arrays[0].clone(), |acc, arr| {
                        acc.iter().zip(arr.iter()).map(|(a, b)| a.min(*b)).collect()
                    })
                }
                AggregationMode::Median => {
                    // Median will be computed per bin after binning
                    // For now, create a placeholder that we'll handle specially
                    arrays[0].clone()
                }
            };

            // Binning
            let chromosome_length = chromosome.length as usize;
            let n_bins = (chromosome_length as f32 / bin_size as f32).ceil() as usize;
            let mut chroms = Vec::with_capacity(n_bins);
            let mut starts = Vec::with_capacity(n_bins);
            let mut ends = Vec::with_capacity(n_bins);
            let mut values: Vec<f32> = Vec::with_capacity(n_bins);

            for bin_idx in 0..n_bins {
                let start = bin_idx * bin_size as usize;
                let end = std::cmp::min((bin_idx + 1) * bin_size as usize, chromosome_length);
                
                // Skip bins that would start at or beyond the chromosome end
                if start >= chromosome_length {
                    continue;
                }
                chroms.push(chromosome.name.clone());
                starts.push(start as u32);
                ends.push(end as u32);

                let bin_value = match aggregation_mode {
                    AggregationMode::Median => {
                        // For median, collect values from all arrays in this bin
                        let mut bin_values: Vec<f32> = Vec::new();
                        for arr in &arrays {
                            if end <= arr.len() {
                                let slice = arr.slice(s![start..end]);
                                bin_values.extend(slice.to_vec());
                            }
                        }
                        if bin_values.is_empty() {
                            0.0
                        } else {
                            bin_values.sort_by(|a, b| {
                                a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                            });
                            let mid = bin_values.len() / 2;
                            if bin_values.len().is_multiple_of(2) {
                                (bin_values[mid - 1] + bin_values[mid]) / 2.0
                            } else {
                                bin_values[mid]
                            }
                        }
                    }
                    _ => {
                        if end <= aggregated.len() {
                            aggregated.slice(s![start..end]).mean().unwrap_or(0.0)
                        } else {
                            0.0
                        }
                    }
                };

                values.push(bin_value);
            }

            let mut df = df![
                "chrom" => chroms,
                "start" => starts,
                "end" => ends,
                "score" => values,
            ]?;
            df.sort_in_place(["chrom", "start"], SortMultipleOptions::default())?;
            Ok(df)
        })
        .collect::<Result<Vec<_>>>()?;

    // Concatenate all DataFrames
    if dfs.is_empty() {
        return Ok(());
    }
    let mut final_df = dfs[0].clone();
    for df in dfs.iter().skip(1) {
        final_df.vstack_mut(df)?;
    }

    // Write chrom sizes to temp file
    let chromsizes_file = tempfile::NamedTempFile::new()?;
    {
        let mut writer = std::io::BufWriter::new(&chromsizes_file);
        for chrom in &chrom_info {
            writeln!(writer, "{}\t{}", chrom.name, chrom.length)?;
        }
    }

    // Write bedgraph to temp file
    info!("Writing temporary bedgraph file");
    let bedgraph_file = tempfile::NamedTempFile::new()?;
    let bedgraph_path = bedgraph_file.path();

    // Use CsvWriter to write tab-separated values
    CsvWriter::new(&bedgraph_file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(&mut final_df)?;

    // Convert to BigWig
    let args = bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs {
        bedgraph: bedgraph_path.to_string_lossy().to_string(),
        chromsizes: chromsizes_file.path().to_string_lossy().to_string(),
        output: output_path.to_string_lossy().to_string(),
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
    info!("Writing output BigWig to {:?}", output_path);
    bigtools::utils::cli::bedgraphtobigwig::bedgraphtobigwig(args)
        .map_err(|e| anyhow::anyhow!("Error converting bedgraph to bigwig: {}", e))?;

    Ok(())
}

#[cfg(test)]

mod tests {
    use super::*;
    use bigtools::BigWigRead;
    use std::collections::HashMap;

    fn create_dummy_bigwig(
        path: &Path,
        chrom_map: HashMap<String, u32>,
        values: Vec<(String, u32, u32, f32)>,
    ) -> Result<()> {
        // Write chrom sizes
        let chromsizes_file = tempfile::NamedTempFile::new()?;
        {
            let mut writer = std::io::BufWriter::new(&chromsizes_file);
            for (chrom, size) in &chrom_map {
                writeln!(writer, "{}\t{}", chrom, size)?;
            }
        }

        // Write bedgraph
        let bedgraph_file = tempfile::NamedTempFile::new()?;
        {
            let mut writer = std::io::BufWriter::new(&bedgraph_file);
            for (chrom, start, end, value) in values {
                writeln!(writer, "{}\t{}\t{}\t{}", chrom, start, end, value)?;
            }
        }

        // Convert
        let args = bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs {
            bedgraph: bedgraph_file.path().to_string_lossy().to_string(),
            chromsizes: chromsizes_file.path().to_string_lossy().to_string(),
            output: path.to_string_lossy().to_string(),
            parallel: "auto".to_string(),
            single_pass: true,
            write_args: bigtools::utils::cli::BBIWriteArgs {
                nthreads: 1,
                nzooms: 0,
                uncompressed: false,
                sorted: "all".to_string(),
                zooms: None,
                block_size: 256,
                items_per_slot: 64,
                inmemory: false,
            },
        };

        bigtools::utils::cli::bedgraphtobigwig::bedgraphtobigwig(args)
            .map_err(|e| anyhow::anyhow!("Error: {}", e))?;

        Ok(())
    }

    #[test]
    fn test_compare_bigwigs_subtraction() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bw1_path = dir.path().join("test1.bw");
        let bw2_path = dir.path().join("test2.bw");
        let out_path = dir.path().join("out.bw");

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 1000);

        let values1 = vec![
            ("chr1".to_string(), 0, 100, 10.0),
            ("chr1".to_string(), 100, 200, 20.0),
        ];
        create_dummy_bigwig(&bw1_path, chrom_map.clone(), values1)?;

        let values2 = vec![
            ("chr1".to_string(), 0, 100, 5.0),
            ("chr1".to_string(), 100, 200, 10.0),
        ];
        create_dummy_bigwig(&bw2_path, chrom_map.clone(), values2)?;

        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_path,
            Comparison::Subtraction,
            10,
            None,
            None,
        )?;

        assert!(out_path.exists());

        // Verify results
        let mut reader = BigWigRead::open_file(out_path.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 200)?
            .collect::<Result<Vec<_>, _>>()?;

        // We expect roughly 5.0 and 10.0 difference
        // Note: binning might affect exact values if not aligned, but here it aligns

        // Check a few points
        let val1 = intervals
            .iter()
            .find(|i| i.start == 0)
            .map(|i| i.value)
            .unwrap();
        assert!((val1 - 5.0).abs() < 1e-5, "Expected 5.0, got {}", val1);

        let val2 = intervals
            .iter()
            .find(|i| i.start == 100)
            .map(|i| i.value)
            .unwrap();
        assert!((val2 - 10.0).abs() < 1e-5, "Expected 10.0, got {}", val2);

        Ok(())
    }

    #[test]
    fn test_compare_bigwigs_ratio() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bw1_path = dir.path().join("test1.bw");
        let bw2_path = dir.path().join("test2.bw");
        let out_path = dir.path().join("out.bw");

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 1000);

        let values1 = vec![("chr1".to_string(), 0, 100, 10.0)];
        create_dummy_bigwig(&bw1_path, chrom_map.clone(), values1)?;

        let values2 = vec![("chr1".to_string(), 0, 100, 5.0)];
        create_dummy_bigwig(&bw2_path, chrom_map.clone(), values2)?;

        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_path,
            Comparison::Ratio,
            10,
            None,
            None,
        )?;

        assert!(out_path.exists());

        // Verify results
        let mut reader = BigWigRead::open_file(out_path.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 100)?
            .collect::<Result<Vec<_>, _>>()?;

        let val1 = intervals
            .iter()
            .find(|i| i.start == 0)
            .map(|i| i.value)
            .unwrap();
        assert!((val1 - 2.0).abs() < 1e-5, "Expected 2.0, got {}", val1);

        Ok(())
    }

    #[test]
    fn test_compare_bigwigs_identical_inputs() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bw1_path = dir.path().join("test1.bw");
        let bw2_path = dir.path().join("test2.bw");
        let out_sub = dir.path().join("out_sub.bw");
        let out_ratio = dir.path().join("out_ratio.bw");
        let out_logratio = dir.path().join("out_logratio.bw");

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 1000);

        let values = vec![
            ("chr1".to_string(), 0, 100, 10.0),
            ("chr1".to_string(), 100, 200, 20.0),
        ];
        create_dummy_bigwig(&bw1_path, chrom_map.clone(), values.clone())?;
        create_dummy_bigwig(&bw2_path, chrom_map.clone(), values)?;

        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_sub,
            Comparison::Subtraction,
            10,
            None,
            None,
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_ratio,
            Comparison::Ratio,
            10,
            None,
            None,
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_logratio,
            Comparison::LogRatio,
            10,
            None,
            None,
        )?;

        // Subtraction should be ~0 everywhere.
        let mut reader = BigWigRead::open_file(out_sub.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 200)?
            .collect::<Result<Vec<_>, _>>()?;
        for iv in &intervals {
            assert!(iv.value.abs() < 1e-5, "Expected ~0, got {}", iv.value);
        }

        // Ratio should be ~1 everywhere.
        let mut reader = BigWigRead::open_file(out_ratio.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 200)?
            .collect::<Result<Vec<_>, _>>()?;
        for iv in &intervals {
            assert!(
                (iv.value - 1.0).abs() < 1e-5,
                "Expected ~1, got {}",
                iv.value
            );
        }

        // LogRatio should be ~0 everywhere.
        let mut reader = BigWigRead::open_file(out_logratio.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 200)?
            .collect::<Result<Vec<_>, _>>()?;
        for iv in &intervals {
            assert!(iv.value.abs() < 1e-5, "Expected ~0, got {}", iv.value);
        }

        Ok(())
    }

    #[test]
    fn test_compare_bigwigs_ratio_and_logratio_are_finite_with_zeros() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bw1_path = dir.path().join("test1.bw");
        let bw2_path = dir.path().join("test2.bw");
        let out_ratio = dir.path().join("out_ratio.bw");
        let out_logratio = dir.path().join("out_logratio.bw");

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 200);

        // Denominator has zeros; we rely on pseudocount to avoid Inf/NaN.
        create_dummy_bigwig(
            &bw1_path,
            chrom_map.clone(),
            vec![("chr1".to_string(), 0, 100, 10.0)],
        )?;
        create_dummy_bigwig(
            &bw2_path,
            chrom_map.clone(),
            vec![("chr1".to_string(), 0, 100, 0.0)],
        )?;

        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_ratio,
            Comparison::Ratio,
            10,
            None,
            Some(1e-3),
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_logratio,
            Comparison::LogRatio,
            10,
            None,
            Some(1e-3),
        )?;

        let mut reader = BigWigRead::open_file(out_ratio.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 100)?
            .collect::<Result<Vec<_>, _>>()?;
        for iv in &intervals {
            assert!(
                iv.value.is_finite(),
                "Ratio should be finite, got {}",
                iv.value
            );
        }

        let mut reader = BigWigRead::open_file(out_logratio.to_str().unwrap())?;
        let intervals: Vec<_> = reader
            .get_interval("chr1", 0, 100)?
            .collect::<Result<Vec<_>, _>>()?;
        for iv in &intervals {
            assert!(
                iv.value.is_finite(),
                "LogRatio should be finite, got {}",
                iv.value
            );
        }

        Ok(())
    }
}
