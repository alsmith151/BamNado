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
///
/// SAFETY: This function is bounds-safe even if `array_chunk.len()` does not match `(end-start)`.
fn read_bp_range_into_array_chunk(
    bw: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    array_chunk: &mut [f32],
) -> Result<()> {
    if array_chunk.is_empty() || end <= start {
        return Ok(());
    }

    // Hard bound: we can only fill as many bp as array_chunk can hold.
    let max_fill_end = start.saturating_add(array_chunk.len() as u32);
    let effective_end = std::cmp::min(end, max_fill_end);
    if effective_end <= start {
        return Ok(());
    }

    let mut bw = bigtools::BigWigRead::open_file(bw)?;
    let values = bw.get_interval(chrom, start, effective_end)?;

    let base_offset = start as usize;
    let chunk_len = array_chunk.len();

    for interval_res in values {
        let interval = interval_res?;

        // Clip interval to requested region
        let s = std::cmp::max(interval.start, start) as usize;
        let e = std::cmp::min(interval.end, effective_end) as usize;
        if s >= e {
            continue;
        }

        for pos in s..e {
            let idx = pos.saturating_sub(base_offset);
            if idx < chunk_len {
                array_chunk[idx] = interval.value;
            }
        }
    }

    Ok(())
}

/// Fill a bp-resolution ndarray array from a BigWig in bounded parallel chunks.
fn fill_array_bp(
    bw_path: &Path,
    chrom: &str,
    chrom_len: usize,
    chunk_size: usize,
    out: &mut Array1<f32>,
) -> Result<()> {
    if out.len() != chrom_len {
        return Err(anyhow::anyhow!(
            "Output array length {} does not match chrom_len {}",
            out.len(),
            chrom_len
        ));
    }

    out.as_slice_mut()
        .expect("contiguous")
        .par_chunks_mut(chunk_size)
        .enumerate()
        .try_for_each(|(chunk_idx, chunk)| -> Result<()> {
            let start = chunk_idx * chunk_size;
            if start >= chrom_len {
                return Ok(());
            }

            // Clamp to both chromosome end and actual chunk length.
            let max_len = chrom_len - start;
            let eff_len = std::cmp::min(chunk.len(), max_len);
            if eff_len == 0 {
                return Ok(());
            }

            let end = start + eff_len;

            read_bp_range_into_array_chunk(
                bw_path,
                chrom,
                start as u32,
                end as u32,
                &mut chunk[..eff_len],
            )?;

            Ok(())
        })?;

    Ok(())
}

/// Fill bp-resolution arrays for many BigWigs in parallel (one array per input file).
fn fill_arrays_bp<P: AsRef<Path> + Send + Sync>(
    bw_paths: &[P],
    chrom: &str,
    chrom_len: usize,
    chunk_size: usize,
) -> Result<Vec<Array1<f32>>> {
    if bw_paths.is_empty() {
        return Ok(vec![]);
    }

    let mut arrays: Vec<Array1<f32>> = (0..bw_paths.len())
        .map(|_| Array1::<f32>::zeros(chrom_len))
        .collect();

    arrays
        .par_iter_mut()
        .zip(bw_paths.par_iter())
        .try_for_each(|(arr, p)| fill_array_bp(p.as_ref(), chrom, chrom_len, chunk_size, arr))?;

    Ok(arrays)
}

/// Turn a single signal into a binned bedGraph-like DataFrame.
fn binned_df_from_signal(chrom: &str, signal: &Array1<f32>, bin_size: u32) -> Result<DataFrame> {
    let len = signal.len();
    if len == 0 {
        return Ok(DataFrame::empty());
    }

    let bin = bin_size as usize;
    if bin == 0 {
        return Err(anyhow::anyhow!("bin_size must be > 0"));
    }

    let n_bins = len.div_ceil(bin);

    let mut chroms = Vec::with_capacity(n_bins);
    let mut starts = Vec::with_capacity(n_bins);
    let mut ends = Vec::with_capacity(n_bins);
    let mut values = Vec::with_capacity(n_bins);

    for i in 0..n_bins {
        let start = i * bin;
        if start >= len {
            break;
        }
        let end = std::cmp::min(start + bin, len);

        chroms.push(chrom.to_string());
        starts.push(start as u32);
        ends.push(end as u32);
        values.push(signal.slice(s![start..end]).mean().unwrap_or(0.0));
    }

    let mut df = df![
        "chrom" => chroms,
        "start" => starts,
        "end" => ends,
        "score" => values,
    ]?;
    df.sort_in_place(["chrom", "start"], SortMultipleOptions::default())?;
    Ok(df)
}

/// Median aggregation per bin across multiple arrays.
/// Semantics: flatten all bp values across all inputs within a bin, then take median.
fn binned_df_median_from_arrays(
    chrom: &str,
    arrays: &[Array1<f32>],
    bin_size: u32,
) -> Result<DataFrame> {
    if arrays.is_empty() {
        return Ok(DataFrame::empty());
    }
    let len = arrays[0].len();
    if len == 0 {
        return Ok(DataFrame::empty());
    }

    let bin = bin_size as usize;
    if bin == 0 {
        return Err(anyhow::anyhow!("bin_size must be > 0"));
    }

    // If lengths disagree, be strict (this should never happen).
    if arrays.iter().any(|a| a.len() != len) {
        return Err(anyhow::anyhow!(
            "Input arrays have different lengths for chromosome {}",
            chrom
        ));
    }

    let n_bins = len.div_ceil(bin);

    let mut chroms = Vec::with_capacity(n_bins);
    let mut starts = Vec::with_capacity(n_bins);
    let mut ends = Vec::with_capacity(n_bins);
    let mut values = Vec::with_capacity(n_bins);

    for i in 0..n_bins {
        let start = i * bin;
        if start >= len {
            break;
        }
        let end = std::cmp::min(start + bin, len);

        chroms.push(chrom.to_string());
        starts.push(start as u32);
        ends.push(end as u32);

        let mut bin_vals = Vec::with_capacity((end - start) * arrays.len());
        for arr in arrays {
            // slice is safe (start/end clamped to len)
            bin_vals.extend(arr.slice(s![start..end]).iter().copied());
        }

        let med = if bin_vals.is_empty() {
            0.0
        } else {
            bin_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mid = bin_vals.len() / 2;
            if bin_vals.len() % 2 == 0 {
                (bin_vals[mid - 1] + bin_vals[mid]) / 2.0
            } else {
                bin_vals[mid]
            }
        };
        values.push(med);
    }

    let mut df = df![
        "chrom" => chroms,
        "start" => starts,
        "end" => ends,
        "score" => values,
    ]?;
    df.sort_in_place(["chrom", "start"], SortMultipleOptions::default())?;
    Ok(df)
}

/// Concatenate per-chromosome DataFrames.
fn concat_dfs(mut dfs: Vec<DataFrame>) -> Result<Option<DataFrame>> {
    if dfs.is_empty() {
        return Ok(None);
    }
    let mut final_df = dfs.remove(0);
    for df in dfs.iter() {
        final_df.vstack_mut(df)?;
    }
    Ok(Some(final_df))
}

/// Write a bedGraph-like DataFrame as a BigWig using bigtools.
fn write_bedgraph_df_to_bigwig(
    df: &mut DataFrame,
    chrom_info: &[bigtools::ChromInfo],
    output_path: &Path,
) -> Result<()> {
    // Chromsizes
    let chromsizes_file = tempfile::NamedTempFile::new()?;
    {
        let mut writer = std::io::BufWriter::new(&chromsizes_file);
        for chrom in chrom_info {
            writeln!(writer, "{}\t{}", chrom.name, chrom.length)?;
        }
    }

    // Bedgraph temp
    info!("Writing temporary bedgraph file");
    let bedgraph_file = tempfile::NamedTempFile::new()?;
    let bedgraph_path = bedgraph_file.path();

    CsvWriter::new(&bedgraph_file)
        .include_header(false)
        .with_separator(b'\t')
        .finish(df)?;

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

/// Compare two BigWig files and output the result as a BigWig.
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
    let chrom_info = bigtools::BigWigRead::open_file(bw1_path.as_ref())?
        .chroms()
        .to_owned()
        .into_iter()
        .sorted_by(|a, b| a.name.cmp(&b.name))
        .collect::<Vec<_>>();

    let chunk_size = chunksize.unwrap_or(1_000_000);
    let pc = pseudocount.unwrap_or(1e-12) as f32;

    let dfs: Vec<DataFrame> = chrom_info
        .iter()
        .progress_with(progress_bar(
            chrom_info.len() as u64,
            "Processing chromosomes".to_string(),
        ))
        .map(|chr| {
            let len = chr.length as usize;

            let mut a1 = Array1::<f32>::zeros(len);
            let mut a2 = Array1::<f32>::zeros(len);

            fill_array_bp(bw1_path.as_ref(), &chr.name, len, chunk_size, &mut a1)?;
            fill_array_bp(bw2_path, &chr.name, len, chunk_size, &mut a2)?;

            let diff = match comparison {
                Comparison::Subtraction => &a1 - &a2,
                Comparison::Ratio => &a1 / (&a2 + pc),
                Comparison::LogRatio => ((&a1 + pc) / (&a2 + pc)).mapv(|x| x.ln()),
            };

            binned_df_from_signal(&chr.name, &diff, bin_size)
        })
        .collect::<Result<Vec<_>>>()?;

    let Some(mut final_df) = concat_dfs(dfs)? else {
        return Ok(());
    };

    write_bedgraph_df_to_bigwig(&mut final_df, &chrom_info, output_path)
}

/// Aggregate multiple BigWigs into a single output BigWig.
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

    let chrom_info = bigtools::BigWigRead::open_file(bw_paths[0].as_ref())?
        .chroms()
        .to_owned()
        .into_iter()
        .sorted_by(|a, b| a.name.cmp(&b.name))
        .collect::<Vec<_>>();

    let chunk_size = chunksize.unwrap_or(1_000_000);
    let pc = pseudocount.unwrap_or(0.0) as f32;

    let dfs: Vec<DataFrame> = chrom_info
        .iter()
        .progress_with(progress_bar(
            chrom_info.len() as u64,
            "Processing chromosomes".to_string(),
        ))
        .map(|chr| {
            let len = chr.length as usize;

            let mut arrays = fill_arrays_bp(bw_paths, &chr.name, len, chunk_size)?;

            // Add pseudocount to all arrays (in-place)
            arrays.iter_mut().for_each(|a| *a = &*a + pc);

            match aggregation_mode {
                AggregationMode::Median => {
                    binned_df_median_from_arrays(&chr.name, &arrays, bin_size)
                }

                AggregationMode::Sum => {
                    let agg = arrays
                        .iter()
                        .skip(1)
                        .fold(arrays[0].clone(), |acc, a| &acc + a);
                    binned_df_from_signal(&chr.name, &agg, bin_size)
                }
                AggregationMode::Mean => {
                    let sum = arrays
                        .iter()
                        .skip(1)
                        .fold(arrays[0].clone(), |acc, a| &acc + a);
                    let agg = sum / arrays.len() as f32;
                    binned_df_from_signal(&chr.name, &agg, bin_size)
                }
                AggregationMode::Max => {
                    let agg = arrays.iter().skip(1).fold(arrays[0].clone(), |acc, a| {
                        acc.iter().zip(a.iter()).map(|(x, y)| x.max(*y)).collect()
                    });
                    binned_df_from_signal(&chr.name, &agg, bin_size)
                }
                AggregationMode::Min => {
                    let agg = arrays.iter().skip(1).fold(arrays[0].clone(), |acc, a| {
                        acc.iter().zip(a.iter()).map(|(x, y)| x.min(*y)).collect()
                    });
                    binned_df_from_signal(&chr.name, &agg, bin_size)
                }
            }
        })
        .collect::<Result<Vec<_>>>()?;

    let Some(mut final_df) = concat_dfs(dfs)? else {
        return Ok(());
    };

    write_bedgraph_df_to_bigwig(&mut final_df, &chrom_info, output_path)
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
