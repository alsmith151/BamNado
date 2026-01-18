use anyhow::Result;
use bigtools::{DEFAULT_BLOCK_SIZE, DEFAULT_ITEMS_PER_SLOT};
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use std::cmp::{Ordering, min};
use std::io::{Read, Seek, Write};
use std::path::Path;

/// Pairwise comparison applied to two binned signals.
#[derive(Debug, Clone, clap::ValueEnum)]
pub enum Comparison {
    /// Signal1 - Signal2.
    Subtraction,
    /// Signal1 / (Signal2 + pseudocount).
    Ratio,
    /// ln((Signal1 + pseudocount) / (Signal2 + pseudocount)).
    LogRatio,
}

/// Reduction applied across multiple BigWig inputs.
#[derive(Debug, Clone, clap::ValueEnum)]
pub enum AggregationMode {
    /// Sum of per-bin means from each input.
    Sum,
    /// Mean of per-bin means from each input.
    Mean,
    /// Per-bin weighted median across all inputs.
    Median,
    /// Per-bin maximum across all inputs.
    Max,
    /// Per-bin minimum across all inputs.
    Min,
}

/// Minimal interval representation (BigWig yields start/end/value).
#[derive(Clone, Copy, Debug)]
struct Iv {
    start: u32,
    end: u32,
    value: f32,
}

/// Load all intervals for a chromosome; ensure sorted by coordinate.
///
/// BigWig intervals are expected to be non-overlapping; this keeps a defensive sort.
fn load_intervals_for_chrom<R: Read + Seek>(
    reader: &mut bigtools::BigWigRead<R>,
    chrom: &str,
    chrom_len: u32,
) -> Result<Vec<Iv>> {
    let it = reader.get_interval(chrom, 0, chrom_len)?;
    let mut out = Vec::new();
    for r in it {
        let iv = r?;
        if iv.end <= iv.start {
            continue;
        }
        out.push(Iv {
            start: iv.start,
            end: iv.end,
            value: iv.value,
        });
    }
    out.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
    Ok(out)
}

/// Compute number of fixed-width bins covering a chromosome.
fn n_bins(chrom_len: u32, bin_size: u32) -> usize {
    let bin = bin_size as usize;
    if bin == 0 {
        return 0;
    }
    (chrom_len as usize).div_ceil(bin)
}

// -----------------------------
// Core abstraction 1: binning
// -----------------------------

/// Accumulates ∫ value(x) dx over fixed-width bins via constant segments,
/// then converts to mean-per-bin by dividing by bin length.
struct BinAccumulator {
    chrom_len: u32,
    bin_size: u32,
    sums: Vec<f64>, // per-bin sum of (value * bp)
}

impl BinAccumulator {
    fn new(chrom_len: u32, bin_size: u32) -> Self {
        Self {
            chrom_len,
            bin_size,
            sums: vec![0.0; n_bins(chrom_len, bin_size)],
        }
    }

    /// Integrate a constant segment [start, end) with given value into bins.
    #[inline]
    fn add_segment(&mut self, start: u32, end: u32, value: f64) {
        if self.bin_size == 0 || end <= start {
            return;
        }

        let start = start.min(self.chrom_len);
        let end = end.min(self.chrom_len);
        if end <= start {
            return;
        }

        let mut pos = start;
        let bin = self.bin_size;

        while pos < end {
            let bi = pos / bin;
            let bin_end = ((bi + 1) * bin).min(self.chrom_len);
            let overlap_end = end.min(bin_end);
            let overlap_len = (overlap_end - pos) as f64;

            // bi is always in-range because sums length == n_bins(chrom_len, bin)
            self.sums[bi as usize] += value * overlap_len;
            pos = overlap_end;
        }
    }

    /// Convert accumulated sums into mean-per-bin values.
    fn into_means(self) -> Vec<f32> {
        let mut out = Vec::with_capacity(self.sums.len());
        for (i, s) in self.sums.into_iter().enumerate() {
            let start = (i as u32) * self.bin_size;
            if start >= self.chrom_len {
                break;
            }
            let end = (start + self.bin_size).min(self.chrom_len);
            let len = (end - start) as f64;
            out.push(if len > 0.0 { (s / len) as f32 } else { 0.0 });
        }
        out
    }
}

// ------------------------------------
// Core abstraction 2: signal stepping
// ------------------------------------

/// Cursor over a piecewise-constant signal defined by non-overlapping intervals with implicit 0 gaps.
///
/// Invariant: at any position `pos`, the signal value is either the active interval's value
/// or 0. The next possible change occurs at either the active interval end or the next interval start.
struct StepSignalCursor<'a> {
    ivs: &'a [Iv],
    idx: usize,
}

impl<'a> StepSignalCursor<'a> {
    fn new(ivs: &'a [Iv]) -> Self {
        Self { ivs, idx: 0 }
    }

    /// Return (value_at_pos, next_change_pos) at `pos`.
    ///
    /// The next change is either the current interval end or the next interval start.
    fn fetch(&mut self, pos: u32, chrom_len: u32) -> (f32, u32) {
        while self.idx < self.ivs.len() && self.ivs[self.idx].end <= pos {
            self.idx += 1;
        }
        if self.idx < self.ivs.len() {
            let iv = self.ivs[self.idx];
            if iv.start <= pos && pos < iv.end {
                (iv.value, iv.end.min(chrom_len))
            } else {
                (0.0, iv.start.min(chrom_len))
            }
        } else {
            (0.0, chrom_len)
        }
    }
}

// -----------------------------
// Binning primitives
// -----------------------------

/// Bin a single BigWig into mean-per-bin (implicit gaps are zero).
fn binned_means_single_bw(
    bw_path: &Path,
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
) -> Result<Vec<f32>> {
    let mut r = bigtools::BigWigRead::open_file(bw_path)?;
    let ivs = load_intervals_for_chrom(&mut r, chrom, chrom_len)?;

    let mut acc = BinAccumulator::new(chrom_len, bin_size);
    for iv in ivs {
        if iv.value == 0.0 {
            continue;
        }
        acc.add_segment(iv.start, iv.end, iv.value as f64);
    }
    Ok(acc.into_means())
}

// -----------------------------
// Compare (2-signal sweep)
// -----------------------------

/// Sweep two signals in lockstep and compute per-bin comparison values.
fn binned_compare_two_bw(
    bw1: &Path,
    bw2: &Path,
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
    comparison: &Comparison,
    pseudocount: f64,
) -> Result<Vec<f32>> {
    let mut r1 = bigtools::BigWigRead::open_file(bw1)?;
    let mut r2 = bigtools::BigWigRead::open_file(bw2)?;
    let a = load_intervals_for_chrom(&mut r1, chrom, chrom_len)?;
    let b = load_intervals_for_chrom(&mut r2, chrom, chrom_len)?;

    let mut c1 = StepSignalCursor::new(&a);
    let mut c2 = StepSignalCursor::new(&b);
    let mut acc = BinAccumulator::new(chrom_len, bin_size);

    let mut pos = 0u32;
    while pos < chrom_len {
        let (v1, n1) = c1.fetch(pos, chrom_len);
        let (v2, n2) = c2.fetch(pos, chrom_len);
        let mut next = n1.min(n2).min(chrom_len);

        // Defensive: ensure forward progress even with malformed inputs.
        if next <= pos {
            next = (pos + 1).min(chrom_len);
        }

        let out_val = match comparison {
            Comparison::Subtraction => (v1 as f64) - (v2 as f64),
            Comparison::Ratio => (v1 as f64) / ((v2 as f64) + pseudocount),
            Comparison::LogRatio => {
                (((v1 as f64) + pseudocount) / ((v2 as f64) + pseudocount)).ln()
            }
        };

        acc.add_segment(pos, next, out_val);
        pos = next;
    }

    Ok(acc.into_means())
}

// -----------------------------
// Extrema (N-signal sweep)
// -----------------------------

/// Sweep N signals and compute a per-bin maximum/minimum.
fn binned_extrema_multi_bw(
    bw_paths: &[&Path],
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
    want_max: bool,
    pseudocount: f64, // added to each signal before extrema (matches original intent)
) -> Result<Vec<f32>> {
    let mut readers = bw_paths
        .iter()
        .map(|p| bigtools::BigWigRead::open_file(*p))
        .collect::<Result<Vec<_>, _>>()?;

    let mut all_ivs: Vec<Vec<Iv>> = Vec::with_capacity(readers.len());
    for r in readers.iter_mut() {
        all_ivs.push(load_intervals_for_chrom(r, chrom, chrom_len)?);
    }

    let mut cursors: Vec<StepSignalCursor<'_>> = all_ivs
        .iter()
        .map(|ivs| StepSignalCursor::new(ivs))
        .collect();
    let mut acc = BinAccumulator::new(chrom_len, bin_size);

    let mut pos = 0u32;
    while pos < chrom_len {
        let mut next = chrom_len;
        let mut ext = if want_max {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };

        for c in cursors.iter_mut() {
            let (v, n) = c.fetch(pos, chrom_len);
            next = next.min(n);
            let vv = (v as f64) + pseudocount;
            ext = if want_max { ext.max(vv) } else { ext.min(vv) };
        }

        if next <= pos {
            next = (pos + 1).min(chrom_len);
        }

        acc.add_segment(pos, next, ext);
        pos = next;
    }

    Ok(acc.into_means())
}

// -----------------------------
// Median (weighted, bin-local)
// -----------------------------

/// Weighted median over a multiset represented as (value, weight=len_in_bp).
fn weighted_median_from_segments(mut pairs: Vec<(f32, u32)>) -> f32 {
    if pairs.is_empty() {
        return 0.0;
    }

    pairs.sort_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let total: u64 = pairs.iter().map(|&(_, w)| w as u64).sum();
    if total == 0 {
        return 0.0;
    }
    let half = total.div_ceil(2); // lower median

    let mut cum: u64 = 0;
    for (val, w) in pairs {
        cum += w as u64;
        if cum >= half {
            return val;
        }
    }
    0.0
}

/// Collect (value, weight) segments covering [bin_start, bin_end) for a single BW.
///
/// Implicit gaps are represented with `gap_value` and weighted by their bp length.
fn collect_weighted_values_for_bw_in_bin(
    ivs: &[Iv],
    bin_start: u32,
    bin_end: u32,
    gap_value: f32,
) -> Vec<(f32, u32)> {
    let mut out: Vec<(f32, u32)> = Vec::new();
    if bin_end <= bin_start {
        return out;
    }

    let mut i = ivs.partition_point(|iv| iv.end <= bin_start);
    let mut pos = bin_start;

    while pos < bin_end {
        if i >= ivs.len() {
            out.push((gap_value, bin_end - pos));
            break;
        }
        let iv = ivs[i];

        // gap before next interval
        if iv.start > pos {
            let gap_end = min(iv.start, bin_end);
            if gap_end > pos {
                out.push((gap_value, gap_end - pos));
                pos = gap_end;
                continue;
            }
        }

        // skip past intervals that end before pos
        if iv.end <= pos {
            i += 1;
            continue;
        }

        // overlap with current interval
        if iv.start <= pos && pos < iv.end {
            let seg_end = min(iv.end, bin_end);
            if seg_end > pos {
                out.push((iv.value, seg_end - pos));
                pos = seg_end;
                if pos >= iv.end {
                    i += 1;
                }
                continue;
            }
        }

        // defensive progress
        if iv.start >= bin_end {
            break;
        }
        i += 1;
    }

    out
}

/// Compute per-bin weighted median across multiple BigWigs.
fn binned_median_multi_bw(
    bw_paths: &[&Path],
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
    pseudocount: f64,
) -> Result<Vec<f32>> {
    let mut readers = bw_paths
        .iter()
        .map(|p| bigtools::BigWigRead::open_file(*p))
        .collect::<Result<Vec<_>, _>>()?;

    let mut all_ivs: Vec<Vec<Iv>> = Vec::with_capacity(readers.len());
    for r in readers.iter_mut() {
        let mut ivs = load_intervals_for_chrom(r, chrom, chrom_len)?;
        if pseudocount != 0.0 {
            for iv in ivs.iter_mut() {
                iv.value = (iv.value as f64 + pseudocount) as f32;
            }
        }
        all_ivs.push(ivs);
    }

    // NOTE: Original per-bp logic added pseudocount to every base, including gaps.
    // To match that exactly, use gap_value = pseudocount as f32.
    let gap_value = pseudocount as f32;

    let nb = n_bins(chrom_len, bin_size);
    let mut out = Vec::with_capacity(nb);

    for bi in 0..nb {
        let start = (bi as u32) * bin_size;
        if start >= chrom_len {
            break;
        }
        let end = min(start + bin_size, chrom_len);

        let mut pairs: Vec<(f32, u32)> = Vec::new();
        for ivs in &all_ivs {
            pairs.extend(collect_weighted_values_for_bw_in_bin(
                ivs, start, end, gap_value,
            ));
        }
        out.push(weighted_median_from_segments(pairs));
    }

    Ok(out)
}

// -----------------------------
// Output writing
// -----------------------------

/// Write bins in bedGraph format, merging consecutive equal-valued bins.
fn write_bedgraph_bins<W: Write>(
    mut w: W,
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
    bins: &[f32],
) -> Result<()> {
    if bin_size == 0 {
        return Err(anyhow::anyhow!("bin_size must be > 0"));
    }

    let mut i = 0usize;
    while i < bins.len() {
        let start = (i as u32) * bin_size;
        if start >= chrom_len {
            break;
        }
        let mut end = min(start + bin_size, chrom_len);
        let v = bins[i];

        // Merge identical consecutive bins to reduce output size.
        let mut j = i + 1;
        while j < bins.len() {
            let s2 = (j as u32) * bin_size;
            if s2 >= chrom_len {
                break;
            }
            let e2 = min(s2 + bin_size, chrom_len);
            if bins[j] != v {
                break;
            }
            end = e2;
            j += 1;
        }

        writeln!(w, "{}\t{}\t{}\t{}", chrom, start, end, v)?;
        i = j;
    }
    Ok(())
}

/// Write per-chromosome bins to a temporary bedGraph and convert to BigWig.
fn write_bedgraph_to_bigwig_from_bins(
    chrom_info: &[bigtools::ChromInfo],
    per_chrom_bins: Vec<(String, Vec<f32>)>,
    bin_size: u32,
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

    // Bedgraph
    info!("Writing temporary bedgraph file");
    let bedgraph_file = tempfile::NamedTempFile::new()?;
    {
        let mut w = std::io::BufWriter::new(&bedgraph_file);
        for (chrom, bins) in per_chrom_bins {
            let clen = chrom_info
                .iter()
                .find(|c| c.name == chrom)
                .map(|c| c.length)
                .unwrap_or(0);
            write_bedgraph_bins(&mut w, &chrom, clen, bin_size, &bins)?;
        }
    }

    // Convert to BigWig
    let args = bigtools::utils::cli::bedgraphtobigwig::BedGraphToBigWigArgs {
        bedgraph: bedgraph_file.path().to_string_lossy().to_string(),
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

// -----------------------------
// Public API
// -----------------------------

/// Compare two BigWigs per-bin and emit a new BigWig.
pub fn compare_bigwigs<P>(
    bw1_path: &P,
    bw2_path: &Path,
    output_path: &Path,
    comparison: Comparison,
    bin_size: u32,
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

    // Small default avoids division-by-zero in ratio/log-ratio modes.
    let pc = pseudocount.unwrap_or(1e-12);

    let mut per_chrom_bins = chrom_info
        .par_iter()
        .enumerate()
        .map(|(idx, chr)| {
            let bins = binned_compare_two_bw(
                bw1_path.as_ref(),
                bw2_path,
                &chr.name,
                chr.length,
                bin_size,
                &comparison,
                pc,
            )?;
            Ok((idx, chr.name.clone(), bins))
        })
        .collect::<Result<Vec<_>>>()?;
    per_chrom_bins.sort_by_key(|(idx, _, _)| *idx);
    let per_chrom_bins = per_chrom_bins
        .into_iter()
        .map(|(_, name, bins)| (name, bins))
        .collect::<Vec<_>>();

    write_bedgraph_to_bigwig_from_bins(&chrom_info, per_chrom_bins, bin_size, output_path)
}

/// Streamed Sum/Mean aggregation: avoids Vec<Vec<_>> and avoids allocating a 2D matrix.
fn binned_sum_or_mean(
    bw_refs: &[&Path],
    chrom: &str,
    chrom_len: u32,
    bin_size: u32,
    pseudocount: f64,
    want_mean: bool,
) -> Result<Vec<f32>> {
    let nb = n_bins(chrom_len, bin_size);
    let mut acc = vec![0f64; nb];

    for p in bw_refs {
        let bins = binned_means_single_bw(p, chrom, chrom_len, bin_size)?;
        for (i, &x) in bins.iter().enumerate() {
            acc[i] += (x as f64) + pseudocount;
        }
    }

    let denom = if want_mean { bw_refs.len() as f64 } else { 1.0 };
    Ok(acc.into_iter().map(|s| (s / denom) as f32).collect())
}

/// Aggregate multiple BigWigs per-bin and emit a new BigWig.
pub fn aggregate_bigwigs<P>(
    bw_paths: &[P],
    output_path: &Path,
    aggregation_mode: AggregationMode,
    bin_size: u32,
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

    let pc = pseudocount.unwrap_or(0.0);
    let bw_refs: Vec<&Path> = bw_paths.iter().map(|p| p.as_ref()).collect();

    let mut per_chrom_bins = chrom_info
        .par_iter()
        .enumerate()
        .map(|(idx, chr)| {
            let clen = chr.length;
            let bins = match aggregation_mode {
                AggregationMode::Sum => {
                    binned_sum_or_mean(&bw_refs, &chr.name, clen, bin_size, pc, false)?
                }
                AggregationMode::Mean => {
                    binned_sum_or_mean(&bw_refs, &chr.name, clen, bin_size, pc, true)?
                }
                AggregationMode::Max => {
                    binned_extrema_multi_bw(&bw_refs, &chr.name, clen, bin_size, true, pc)?
                }
                AggregationMode::Min => {
                    binned_extrema_multi_bw(&bw_refs, &chr.name, clen, bin_size, false, pc)?
                }
                AggregationMode::Median => {
                    binned_median_multi_bw(&bw_refs, &chr.name, clen, bin_size, pc)?
                }
            };
            Ok((idx, chr.name.clone(), bins))
        })
        .collect::<Result<Vec<_>>>()?;
    per_chrom_bins.sort_by_key(|(idx, _, _)| *idx);
    let per_chrom_bins = per_chrom_bins
        .into_iter()
        .map(|(_, name, bins)| (name, bins))
        .collect::<Vec<_>>();

    write_bedgraph_to_bigwig_from_bins(&chrom_info, per_chrom_bins, bin_size, output_path)
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
    fn test_write_bedgraph_bins_collapses_equal_signals() -> Result<()> {
        let tmp = tempfile::NamedTempFile::new()?;
        {
            let mut w = std::io::BufWriter::new(&tmp);
            let bins = vec![1.0f32, 1.0, 2.0, 2.0, 1.0];
            write_bedgraph_bins(&mut w, "chr1", 500, 100, &bins)?;
        }
        let contents = std::fs::read_to_string(tmp.path())?;
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines.len(), 3);

        let parse_line = |s: &str| -> Result<(String, u32, u32, f32)> {
            let parts: Vec<&str> = s.split_whitespace().collect();
            Ok((
                parts[0].to_string(),
                parts[1].parse()?,
                parts[2].parse()?,
                parts[3].parse()?,
            ))
        };

        let l0 = parse_line(lines[0])?;
        assert_eq!(l0.0, "chr1");
        assert_eq!(l0.1, 0);
        assert_eq!(l0.2, 200);
        assert!((l0.3 - 1.0).abs() < 1e-6);

        let l1 = parse_line(lines[1])?;
        assert_eq!(l1.1, 200);
        assert_eq!(l1.2, 400);
        assert!((l1.3 - 2.0).abs() < 1e-6);

        let l2 = parse_line(lines[2])?;
        assert_eq!(l2.1, 400);
        assert_eq!(l2.2, 500);
        assert!((l2.3 - 1.0).abs() < 1e-6);

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

        compare_bigwigs(&bw1_path, &bw2_path, &out_path, Comparison::Ratio, 10, None)?;

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
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_ratio,
            Comparison::Ratio,
            10,
            None,
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_logratio,
            Comparison::LogRatio,
            10,
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
            Some(1e-3),
        )?;
        compare_bigwigs(
            &bw1_path,
            &bw2_path,
            &out_logratio,
            Comparison::LogRatio,
            10,
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
