use anyhow::{Result, bail};
use std::collections::HashSet;
use std::io::{Read, Seek};
use std::path::Path;

/// Inferred normalisation method.
#[derive(Debug, Clone, PartialEq)]
pub enum NormMethod {
    Cpm,
    Rpkm,
    /// Both CPM and RPKM give a plausible library size — cannot distinguish.
    Ambiguous,
    /// Neither CPM nor RPKM gives a plausible library size (RPGC, BPM, raw, …).
    Unknown,
}

impl std::fmt::Display for NormMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NormMethod::Cpm => write!(f, "CPM"),
            NormMethod::Rpkm => write!(f, "RPKM"),
            NormMethod::Ambiguous => write!(f, "Ambiguous (CPM or RPKM)"),
            NormMethod::Unknown => write!(f, "Unknown"),
        }
    }
}

#[derive(Debug)]
pub struct Warnings {
    /// More than 2 distinct interval widths — variable bin sizes found.
    pub variable_bin_sizes: bool,
    /// second_min / min is not near an integer — smoothing / kernel applied.
    pub smoothing_detected: bool,
    /// ratio > 2 and near integer — pseudocount shifted the signal floor.
    pub pseudocount_detected: bool,
    /// Minimum value came from a chromosome that is < 0.1 % of genome size.
    pub min_from_small_chrom: bool,
}

#[derive(Debug)]
pub struct InferScaleResult {
    /// Inferred normalisation method.
    pub norm_method: NormMethod,
    /// Scale factor s such that raw ≈ normalised × s.
    pub scale_factor: f64,
    /// Implied library size (total mapped reads or read-pairs).
    pub library_size: f64,
    /// Canonical bin size (bp) inferred from interval widths.
    pub bin_size: u32,
    /// Global minimum non-zero bin value.
    pub min_val: f64,
    /// Second-smallest distinct non-zero bin value.
    pub second_min_val: f64,
    /// second_min / min (diagnostic ratio).
    pub ratio: f64,
    /// Estimated pseudocount (if detected).
    pub pseudocount: Option<f64>,
    /// True when the confidence criterion was met before exhausting chromosomes.
    pub confident: bool,
    pub warnings: Warnings,
    /// Number of chromosomes scanned.
    pub chroms_scanned: usize,
}

/// Infer scale factor and normalisation method from a BigWig file.
pub fn infer_scale_factor(bw_path: &Path) -> Result<InferScaleResult> {
    let mut reader = bigtools::BigWigRead::open_file(bw_path)?;
    let chrom_info: Vec<bigtools::ChromInfo> = reader.chroms().to_owned();

    if chrom_info.is_empty() {
        bail!("BigWig file has no chromosomes");
    }

    let genome_size: u64 = chrom_info.iter().map(|c| c.length as u64).sum();

    // Sort: mid-sized autosomes first, large deferred, small/unplaced last.
    let mut sorted_chroms = chrom_info.clone();
    sorted_chroms.sort_by_key(|c| chrom_priority(&c.name, c.length, genome_size));

    let mut global_min = f64::MAX;
    let mut global_second_min = f64::MAX;
    let mut global_starts: HashSet<u32> = HashSet::new();
    let mut min_confirmations: u32 = 0;
    let mut chroms_scanned: usize = 0;
    let mut confident = false;
    let mut min_chrom_name = String::new();

    for chrom in &sorted_chroms {
        let (cmin, c2min, cwidths) = scan_chrom(&mut reader, &chrom.name, chrom.length)?;
        chroms_scanned += 1;
        global_starts.extend(cwidths);

        if cmin == f64::MAX {
            continue;
        }

        const EPS: f64 = 1e-9;

        if cmin < global_min - EPS {
            let old_min = global_min;
            let old_second = global_second_min;
            // Candidates for new second_min: old global_min, old second_min, this chrom's second_min.
            global_second_min = [old_min, old_second, c2min]
                .iter()
                .copied()
                .filter(|&v| v > cmin + EPS)
                .fold(f64::MAX, f64::min);
            global_min = cmin;
            min_chrom_name = chrom.name.clone();
            min_confirmations = 1;
        } else if (cmin - global_min).abs() < EPS {
            min_confirmations += 1;
            if c2min < global_second_min - EPS {
                global_second_min = c2min;
            }
        } else {
            // cmin > global_min: update second_min only.
            if cmin < global_second_min - EPS {
                global_second_min = cmin;
            }
            if c2min > global_min + EPS && c2min < global_second_min - EPS {
                global_second_min = c2min;
            }
        }

        if min_confirmations >= 2 && global_second_min < f64::MAX {
            let ratio = global_second_min / global_min;
            if (ratio - ratio.round()).abs() < 0.01 {
                confident = true;
                break;
            }
        }
    }

    if global_min == f64::MAX {
        bail!("No non-zero values found in BigWig file");
    }

    // GCD of all observed widths (end-truncated bins are excluded in scan_chrom).
    // GCD = 0 → no widths collected; GCD = 1 → no common factor → truly mixed bin sizes.
    let canonical_bin_size = global_starts.iter().copied().fold(0u32, gcd);
    let variable_bin_sizes = canonical_bin_size <= 1;

    let ratio = if global_second_min < f64::MAX {
        global_second_min / global_min
    } else {
        f64::NAN
    };

    let smoothing_detected = !ratio.is_nan() && (ratio - ratio.round()).abs() >= 0.01;

    // Pseudocount: ratio > 2 and near integer.
    // Min bin = pseudocount-only bin (0 raw reads). s_corrected = 1 / second_min.
    let (pseudocount, scale_factor) =
        if !ratio.is_nan() && ratio > 2.0 + 0.05 && (ratio - ratio.round()).abs() < 0.05 {
            let p = 1.0 / (ratio - 1.0);
            (Some(p), 1.0 / global_second_min)
        } else {
            (None, 1.0 / global_min)
        };

    // CPM: min_val = 1e6 / N  →  N = 1e6 / min_val
    // RPKM: min_val = 1e9 / (N × bin_size)  →  N = 1e9 / (min_val × bin_size)
    let n_cpm = 1e6 / global_min;
    let n_rpkm = if canonical_bin_size > 0 {
        1e9 / (global_min * canonical_bin_size as f64)
    } else {
        f64::NAN
    };

    let plausible = |n: f64| (1e6..=2e9).contains(&n);
    let norm_method = match (plausible(n_cpm), !n_rpkm.is_nan() && plausible(n_rpkm)) {
        (true, false) => NormMethod::Cpm,
        (false, true) => NormMethod::Rpkm,
        (true, true) => NormMethod::Ambiguous,
        (false, false) => NormMethod::Unknown,
    };

    let library_size = match norm_method {
        NormMethod::Rpkm => n_rpkm,
        _ => n_cpm,
    };

    let min_chrom_fraction = sorted_chroms
        .iter()
        .find(|c| c.name == min_chrom_name)
        .map(|c| c.length as f64 / genome_size as f64)
        .unwrap_or(1.0);
    let min_from_small_chrom = min_chrom_fraction < 0.001;

    Ok(InferScaleResult {
        norm_method,
        scale_factor,
        library_size,
        bin_size: canonical_bin_size,
        min_val: global_min,
        second_min_val: global_second_min,
        ratio,
        pseudocount,
        confident,
        warnings: Warnings {
            variable_bin_sizes,
            smoothing_detected,
            pseudocount_detected: pseudocount.is_some(),
            min_from_small_chrom,
        },
        chroms_scanned,
    })
}

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn scan_chrom<R: Read + Seek>(
    reader: &mut bigtools::BigWigRead<R>,
    chrom: &str,
    chrom_len: u32,
) -> Result<(f64, f64, HashSet<u32>)> {
    let it = reader.get_interval(chrom, 0, chrom_len)?;
    let mut min1 = f64::MAX;
    let mut min2 = f64::MAX;
    let mut starts: HashSet<u32> = HashSet::new();

    for r in it {
        let iv = r?;
        if iv.value <= 0.0 || iv.end <= iv.start {
            continue;
        }
        let v = iv.value as f64;

        // Track non-zero start positions for bin_size inference.
        // All intervals in a fixed-bin BigWig start at multiples of bin_size, including
        // end-truncated ones — so starts are robust to chromosome boundary truncation.
        // start=0 contributes nothing to GCD (identity), so skip it.
        if iv.start > 0 {
            starts.insert(iv.start);
        }

        if v < min1 {
            if min1 < min2 {
                min2 = min1;
            }
            min1 = v;
        } else if v > min1 + 1e-9 && v < min2 {
            min2 = v;
        }
    }

    Ok((min1, min2, starts))
}

/// Returns sort key: 0 = mid-sized autosomes, 1 = large, 2 = small/unplaced/mito.
pub(crate) fn chrom_priority(name: &str, len: u32, genome_size: u64) -> u8 {
    if genome_size == 0 {
        return 1;
    }
    let fraction = len as f64 / genome_size as f64;
    let lower = name.to_ascii_lowercase();

    if lower == "chrm"
        || lower == "mt"
        || lower == "chrmt"
        || lower == "m"
        || lower.contains('_')
        || lower.contains("random")
        || lower.contains("un")
        || fraction < 0.001
    {
        return 2;
    }

    if (0.02..=0.25).contains(&fraction) {
        0
    } else if fraction > 0.25 {
        1
    } else {
        2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bigtools::beddata::BedParserStreamingIterator;
    use bigtools::{BigWigWrite, Value};
    use std::collections::HashMap;

    // Standard test parameters used throughout: bin_size=10, ~10M reads.
    // CPM min_val = 1e6/N.  RPKM min_val = 1e9/(N*bin_size).

    fn make_bigwig(
        path: &std::path::Path,
        chroms: HashMap<String, u32>,
        intervals: Vec<(String, u32, u32, f32)>,
    ) -> anyhow::Result<()> {
        let iter = intervals.iter().map(|(c, s, e, v)| {
            (
                c.as_str(),
                Value {
                    start: *s,
                    end: *e,
                    value: *v,
                },
            )
        });
        let bed_iter = BedParserStreamingIterator::wrap_infallible_iter(iter, true);
        let writer = BigWigWrite::create_file(path, chroms)?;
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()?;
        writer.write(bed_iter, runtime)?;
        Ok(())
    }

    // ─── chrom_priority unit tests ───────────────────────────────────────────

    #[test]
    fn test_priority_mid_sized_autosome() {
        // 10% of genome → mid-sized
        assert_eq!(chrom_priority("chr5", 1000, 10000), 0);
    }

    #[test]
    fn test_priority_large_chrom() {
        // 30% of genome → large
        assert_eq!(chrom_priority("chr1", 3000, 10000), 1);
    }

    #[test]
    fn test_priority_small_chrom() {
        // 1% of genome, no special name → small
        assert_eq!(chrom_priority("chr22", 100, 10000), 2);
    }

    #[test]
    fn test_priority_mito_name() {
        // chrM is always deprioritised regardless of fraction
        assert_eq!(chrom_priority("chrM", 1000, 10000), 2);
        assert_eq!(chrom_priority("MT", 1000, 10000), 2);
    }

    #[test]
    fn test_priority_unplaced_scaffold() {
        // Underscore in name → unplaced
        assert_eq!(chrom_priority("chr1_random", 1000, 10000), 2);
        assert_eq!(chrom_priority("chrUn_gl000220", 500, 10000), 2);
    }

    #[test]
    fn test_priority_zero_genome_size() {
        // Degenerate — return 1 (not 0, not 2)
        assert_eq!(chrom_priority("chr1", 1000, 0), 1);
    }

    #[test]
    fn test_priority_tiny_fraction() {
        // < 0.1% of genome → deprioritised even without special name
        assert_eq!(chrom_priority("chrSmall", 1, 1_000_000), 2);
    }

    // ─── infer_scale_factor integration tests ────────────────────────────────

    /// Build a BigWig where the same minimum value appears on two chromosomes.
    /// Standard parameters: bin_size=10, N_total=10M reads.
    fn make_two_chrom_bigwig(
        dir: &tempfile::TempDir,
        min_val: f32,
        second_val: f32,
        bin_size: u32,
    ) -> anyhow::Result<std::path::PathBuf> {
        let path = dir.path().join("test.bw");
        let chrom_len: u32 = bin_size * 20; // 20 bins each
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), chrom_len);
        chroms.insert("chr7".to_string(), chrom_len);

        let mut intervals: Vec<(String, u32, u32, f32)> = Vec::new();
        for chrom in &["chr5", "chr7"] {
            // One singleton bin (min_val), rest doublet (second_val)
            intervals.push((chrom.to_string(), 0, bin_size, min_val));
            for i in 1..20u32 {
                intervals.push((
                    chrom.to_string(),
                    i * bin_size,
                    (i + 1) * bin_size,
                    second_val,
                ));
            }
        }
        // Intervals must be sorted per chrom — they already are.
        make_bigwig(&path, chroms, intervals)?;
        Ok(path)
    }

    #[test]
    fn test_cpm_detection() -> anyhow::Result<()> {
        // N=100M, bin_size=10: min_val=0.01, n_cpm=100M ✓, n_rpkm=10B ✗ → CPM
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path)?;

        assert_eq!(r.norm_method, NormMethod::Cpm);
        assert!((r.scale_factor - 100.0).abs() < 1.0, "scale_factor ≈ 100");
        assert!((r.library_size - 100_000_000.0).abs() < 1_000_000.0);
        assert_eq!(r.bin_size, 10);
        assert!((r.ratio - 2.0).abs() < 0.01);
        assert!(r.confident);
        assert!(r.pseudocount.is_none());
        assert!(!r.warnings.smoothing_detected);
        assert!(!r.warnings.pseudocount_detected);
        Ok(())
    }

    #[test]
    fn test_rpkm_detection() -> anyhow::Result<()> {
        // N=10M, bin_size=10: min_val=10.0, n_cpm=100K ✗, n_rpkm=10M ✓ → RPKM
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 10.0, 20.0, 10)?;
        let r = infer_scale_factor(&path)?;

        assert_eq!(r.norm_method, NormMethod::Rpkm);
        assert!((r.scale_factor - 0.1).abs() < 0.001, "scale_factor ≈ 0.1");
        assert!((r.library_size - 10_000_000.0).abs() < 100_000.0);
        assert_eq!(r.bin_size, 10);
        assert!((r.ratio - 2.0).abs() < 0.01);
        assert!(r.confident);
        Ok(())
    }

    #[test]
    fn test_raw_unknown_detection_exits_quickly() -> anyhow::Result<()> {
        // Raw count data: min_val=1000.0 (not plausible for CPM or RPKM).
        // Confident after 2 chroms → early exit without scanning remaining genome.
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 1000.0, 2000.0, 10)?;
        let r = infer_scale_factor(&path)?;

        assert_eq!(r.norm_method, NormMethod::Unknown);
        assert!(r.confident, "should exit confidently after 2 chroms");
        assert_eq!(
            r.chroms_scanned, 2,
            "early exit after 2 confirmation chroms"
        );
        assert!((r.ratio - 2.0).abs() < 0.01);
        Ok(())
    }

    #[test]
    fn test_pseudocount_detected() -> anyhow::Result<()> {
        // p=0.5: ratio = (1+p)/p = 3.0 > 2 and integer.
        // min = p/N*1e6 = 0.5/8M*1e6 = 0.0625 (exact f32 = 1/16)
        // second = (1+p)/N*1e6 = 1.5/8M*1e6 = 0.1875 (exact f32 = 3/16)
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.0625, 0.1875, 10)?;
        let r = infer_scale_factor(&path)?;

        assert!(
            r.warnings.pseudocount_detected,
            "pseudocount should be flagged"
        );
        let p = r.pseudocount.expect("pseudocount should be Some");
        assert!((p - 0.5).abs() < 0.01, "p ≈ 0.5, got {p}");
        // Corrected scale_factor = 1/second_min = 1/0.1875 ≈ 5.333
        assert!(
            (r.scale_factor - (1.0 / 0.1875)).abs() < 0.01,
            "s_corrected = 1/second_min"
        );
        assert!((r.ratio - 3.0).abs() < 0.01);
        assert!(r.confident);
        Ok(())
    }

    #[test]
    fn test_smoothing_detected() -> anyhow::Result<()> {
        // Non-integer ratio (1.5): fractional read contributions / kernel smoothing.
        // min=0.25 (exact f32), second=0.375 (exact f32=3/8), ratio=1.5
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.25, 0.375, 10)?;
        let r = infer_scale_factor(&path)?;

        assert!(r.warnings.smoothing_detected, "smoothing should be flagged");
        assert!(!r.confident, "non-integer ratio → not confident");
        assert!((r.ratio - 1.5).abs() < 0.05);
        Ok(())
    }

    #[test]
    fn test_variable_bin_sizes_warning() -> anyhow::Result<()> {
        // Non-zero start positions 10 and 21 → GCD(10,21)=1 → variable_bin_sizes warning.
        // Simulates a file where two incommensurable bin grids are mixed.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("varbins.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 50u32);
        chroms.insert("chr7".to_string(), 50u32);
        // Non-zero starts visible after skipping start=0: {10, 21} → GCD=1 → variable.
        let intervals = vec![
            ("chr5".to_string(), 0u32, 10u32, 0.01f32),
            ("chr5".to_string(), 10, 21, 0.02),
            ("chr5".to_string(), 21, 50, 0.02),
            ("chr7".to_string(), 0, 10, 0.01),
            ("chr7".to_string(), 10, 21, 0.02),
            ("chr7".to_string(), 21, 50, 0.02),
        ];
        make_bigwig(&path, chroms, intervals)?;
        let r = infer_scale_factor(&path)?;

        assert!(
            r.warnings.variable_bin_sizes,
            "coprime starts should flag variable_bin_sizes"
        );
        Ok(())
    }

    #[test]
    fn test_collapsed_bins_not_variable() -> anyhow::Result<()> {
        // Collapsed BigWig: widths 10, 20, 30 (all ×10), starts {10, 30, 60}.
        // GCD({10,30,60}) = 10 → canonical_bin_size=10, variable_bin_sizes=false.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("collapsed.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 200u32);
        chroms.insert("chr7".to_string(), 200u32);
        let mut intervals = Vec::new();
        for chrom in &["chr5", "chr7"] {
            intervals.push((chrom.to_string(), 0u32, 10u32, 0.01f32));
            intervals.push((chrom.to_string(), 10, 30, 0.02));
            intervals.push((chrom.to_string(), 30, 60, 0.03));
            intervals.push((chrom.to_string(), 60, 200, 0.02));
        }
        make_bigwig(&path, chroms, intervals)?;
        let r = infer_scale_factor(&path)?;

        assert!(
            !r.warnings.variable_bin_sizes,
            "all-multiple starts should not trigger variable warning"
        );
        assert_eq!(r.bin_size, 10, "GCD of starts {{10,30,60}} = 10");
        Ok(())
    }

    #[test]
    fn test_confident_after_two_chroms() -> anyhow::Result<()> {
        // Same min on 2 chroms + integer ratio → confident=true, scan stops at 2.
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path)?;

        assert!(r.confident);
        assert_eq!(r.chroms_scanned, 2);
        Ok(())
    }

    #[test]
    fn test_not_confident_single_chrom() -> anyhow::Result<()> {
        // Only one chromosome → confirmations=1 → confident=false.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("single.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 200u32);
        let intervals = vec![
            ("chr5".to_string(), 0u32, 10u32, 0.01f32),
            ("chr5".to_string(), 10, 20, 0.02),
        ];
        make_bigwig(&path, chroms, intervals)?;
        let r = infer_scale_factor(&path)?;

        assert!(!r.confident);
        assert_eq!(r.chroms_scanned, 1);
        Ok(())
    }

    #[test]
    fn test_scale_factor_round_trip_cpm() {
        // Verify the math: raw = cpm_val * scale_factor recovers original count.
        // N=100M, 1 read in a 10 bp bin → CPM = 1e6/100M = 0.01
        let cpm_val: f64 = 0.01;
        let scale = 1.0 / cpm_val; // = 100.0 = N/1e6
        let raw = cpm_val * scale;
        assert!((raw - 1.0).abs() < 1e-10); // recovers exactly 1 read
    }

    #[test]
    fn test_scale_factor_round_trip_rpkm() {
        // N=10M, 1 read in a 10 bp bin → RPKM = 1e9/(N*bin) = 1e9/(10M*10) = 10.0
        let rpkm_val: f64 = 10.0;
        let scale = 1.0 / rpkm_val; // = 0.1 = N*bin/1e9
        // raw_read_density = rpkm_val * scale = 1.0  (1 read per bin)
        let raw = rpkm_val * scale;
        assert!((raw - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_on_empty_bigwig() {
        // BigWig writer refuses an empty interval list → make_bigwig itself errors.
        // infer_scale_factor would also error ("no non-zero values") if given such a file.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 1000u32);
        let result = make_bigwig(&path, chroms, vec![]);
        assert!(result.is_err(), "bigtools rejects empty interval list");
    }
}
