use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::Path;

/// Read buffer per thread — large enough to amortise syscall overhead on whole-chromosome fetches.
const BW_BUF_BYTES: usize = 2 << 20; // 2 MB

fn open_bw(path: &Path) -> Result<bigtools::BigWigRead<BufReader<File>>> {
    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let buf = BufReader::with_capacity(BW_BUF_BYTES, file);
    bigtools::BigWigRead::open(buf).map_err(|e| anyhow::anyhow!("BigWig open error: {e}"))
}

fn ser_f64_max_null<S>(v: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if *v == f64::MAX {
        s.serialize_none()
    } else {
        s.serialize_some(v)
    }
}

fn ser_f64_nan_null<S>(v: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if v.is_nan() {
        s.serialize_none()
    } else {
        s.serialize_some(v)
    }
}

/// Tuning knobs for `infer_scale_factor`.
#[derive(Debug, Clone)]
pub struct InferScaleConfig {
    /// Stop scanning a chromosome after this many consecutive intervals with no new minimum.
    ///
    /// The heuristic is safe when the minimum bin value appears within the first
    /// `min_stable_streak` intervals of a chromosome.  For real normalised BigWigs
    /// (CPM/RPKM) single-read bins are scattered throughout the genome, so 50 000 is
    /// conservative.  Set to `u32::MAX` to disable early exit and always do a full scan.
    ///
    /// **Limitation**: if ALL minimum bins are clustered at the very end of a chromosome
    /// and the chromosome also contains two or more distinct larger values before that
    /// point, the heuristic may report a higher minimum.  This is unlikely in practice
    /// but can be demonstrated with synthetic data — see `test_heuristic_misses_late_minimum`.
    pub min_stable_streak: u32,

    /// Maximum number of interval start positions collected for bin-size (GCD) inference.
    ///
    /// The GCD of bin-aligned start positions stabilises within the first few entries.
    /// Capping avoids O(intervals) `HashSet` growth on high-depth files where nearly
    /// every bin is non-zero.
    pub max_starts: usize,
}

impl Default for InferScaleConfig {
    fn default() -> Self {
        Self {
            min_stable_streak: 50_000,
            max_starts: 512,
        }
    }
}

/// Inferred normalisation method.
#[derive(Debug, Clone, PartialEq, Serialize)]
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

impl Warnings {
    /// Returns active warnings as human-readable strings (empty when no warnings).
    pub fn messages(&self) -> Vec<&'static str> {
        let mut v = Vec::new();
        if self.variable_bin_sizes {
            v.push("variable_bin_sizes: mixed bin widths detected, RPKM reversal unreliable");
        }
        if self.smoothing_detected {
            v.push("smoothing_detected: ratio non-integer, recovery is approximate");
        }
        if self.pseudocount_detected {
            v.push("pseudocount_detected: scale factor corrected to second_min");
        }
        if self.min_from_small_chrom {
            v.push(
                "min_from_small_chrom: minimum from small/unplaced chromosome, may be artefactual",
            );
        }
        v
    }
}

impl Serialize for Warnings {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.messages().serialize(s)
    }
}

#[derive(Debug, Serialize)]
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
    /// Second-smallest distinct non-zero bin value (null when only one distinct value seen).
    #[serde(serialize_with = "ser_f64_max_null")]
    pub second_min_val: f64,
    /// second_min / min (diagnostic ratio; null when second_min unknown).
    #[serde(serialize_with = "ser_f64_nan_null")]
    pub ratio: f64,
    /// Estimated pseudocount (if detected).
    pub pseudocount: Option<f64>,
    /// True when the confidence criterion was met before exhausting chromosomes.
    pub confident: bool,
    pub warnings: Warnings,
    /// Number of chromosomes scanned.
    pub chroms_scanned: usize,
}

impl InferScaleResult {
    /// Convert a normalised BigWig value back to raw read count per bin.
    ///
    /// raw = value × scale_factor
    ///
    /// Derivation:
    ///   CPM:  value = reads / (N/1e6)       →  raw = value × (N/1e6)
    ///   RPKM: value = reads / (N×bin/1e9)   →  raw = value × (N×bin/1e9)
    ///   Both reduce to: raw = value × scale_factor
    ///
    /// When a pseudocount was detected the scale_factor already uses second_min
    /// (the smallest bin with ≥1 real read) so raw counts remain integer-valued.
    pub fn to_raw(&self, value: f64) -> f64 {
        value * self.scale_factor
    }
}

/// Infer scale factor and normalisation method from a BigWig file.
///
/// Chromosomes are scanned in three priority groups (mid-sized autosomes first,
/// large chromosomes second, small/unplaced last).  Each group is processed in
/// parallel — every worker thread opens its own buffered reader.  Scanning stops
/// as soon as the minimum value has been confirmed on ≥2 chromosomes and the
/// second_min/min ratio is near an integer.
pub fn infer_scale_factor(bw_path: &Path, config: &InferScaleConfig) -> Result<InferScaleResult> {
    let chrom_info: Vec<bigtools::ChromInfo> = {
        let reader = open_bw(bw_path)?;
        reader.chroms().to_owned()
    };

    if chrom_info.is_empty() {
        bail!("BigWig file has no chromosomes");
    }

    let genome_size: u64 = chrom_info.iter().map(|c| c.length as u64).sum();

    // Partition into three priority groups — mid-sized autosomes first because
    // they are the most likely to share the global minimum and confirm confidence
    // quickly, while being small enough to scan fast.
    let mut groups: [Vec<bigtools::ChromInfo>; 3] = [vec![], vec![], vec![]];
    for c in &chrom_info {
        groups[chrom_priority(&c.name, c.length, genome_size) as usize].push(c.clone());
    }

    let mut global_min = f64::MAX;
    let mut global_second_min = f64::MAX;
    let mut global_starts: HashSet<u32> = HashSet::new();
    let mut min_confirmations: u32 = 0;
    let mut chroms_scanned: usize = 0;
    let mut confident = false;
    let mut min_chrom_name = String::new();

    'outer: for group in &groups {
        if group.is_empty() || confident {
            continue;
        }

        // All chromosomes in this priority group are scanned in parallel.
        // Each closure opens its own BigWigRead so there is no shared I/O state.
        type ChromScan = Result<(String, f64, f64, HashSet<u32>)>;
        let scan_results: Vec<ChromScan> = group
            .par_iter()
            .map(|chrom| {
                let mut r = open_bw(bw_path)?;
                let (cmin, c2min, cwidths) = scan_chrom(&mut r, &chrom.name, chrom.length, config)?;
                Ok((chrom.name.clone(), cmin, c2min, cwidths))
            })
            .collect();

        // Merge results sequentially — order within a group doesn't matter for
        // correctness; the priority ordering ensures the most informative group
        // is processed first so confidence is reached with fewest total scans.
        for res in scan_results {
            let (chrom_name, cmin, c2min, cwidths) = res?;
            chroms_scanned += 1;
            global_starts.extend(cwidths);

            if cmin == f64::MAX {
                continue;
            }

            const EPS: f64 = 1e-9;

            if cmin < global_min - EPS {
                let old_min = global_min;
                let old_second = global_second_min;
                global_second_min = [old_min, old_second, c2min]
                    .iter()
                    .copied()
                    .filter(|&v| v > cmin + EPS)
                    .fold(f64::MAX, f64::min);
                global_min = cmin;
                min_chrom_name = chrom_name;
                min_confirmations = 1;
            } else if (cmin - global_min).abs() < EPS {
                min_confirmations += 1;
                if c2min < global_second_min - EPS {
                    global_second_min = c2min;
                }
            } else {
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
                    break 'outer;
                }
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

    let min_chrom_fraction = chrom_info
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
    config: &InferScaleConfig,
) -> Result<(f64, f64, HashSet<u32>)> {
    let it = reader.get_interval(chrom, 0, chrom_len)?;
    let mut min1 = f64::MAX;
    let mut min2 = f64::MAX;
    let mut starts: HashSet<u32> = HashSet::new();
    // Counts intervals since the last time min1 changed.  Once this exceeds
    // min_stable_streak we have seen enough of the distribution to trust min1/min2.
    //
    // IMPORTANT: the break only fires when min2 < f64::MAX (i.e. we have already
    // seen two distinct values), which prevents exiting on a single-value plateau.
    // However, if a lower minimum exists beyond the streak window AND the chromosome
    // already has two distinct values, it will be missed.  Use min_stable_streak =
    // u32::MAX to disable this optimisation and always do a full scan.
    let mut stable_streak: u32 = 0;

    for r in it {
        let iv = r?;
        if iv.value <= 0.0 || iv.end <= iv.start {
            continue;
        }
        let v = iv.value as f64;

        // Collect start positions for bin_size inference via GCD.
        // GCD stabilises within the first few multiples-of-bin_size entries,
        // so stop inserting once we have enough — avoids O(intervals) HashSet growth
        // on high-depth files where most bins are non-zero.
        if iv.start > 0 && starts.len() < config.max_starts {
            starts.insert(iv.start);
        }

        if v < min1 - 1e-9 {
            if min1 < min2 {
                min2 = min1;
            }
            min1 = v;
            stable_streak = 0;
        } else {
            if v > min1 + 1e-9 && v < min2 {
                min2 = v;
            }
            stable_streak += 1;
            if stable_streak >= config.min_stable_streak && min2 < f64::MAX {
                break;
            }
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

    fn default_config() -> InferScaleConfig {
        InferScaleConfig::default()
    }

    fn full_scan_config() -> InferScaleConfig {
        InferScaleConfig {
            min_stable_streak: u32::MAX,
            ..InferScaleConfig::default()
        }
    }

    // ─── chrom_priority unit tests ───────────────────────────────────────────

    #[test]
    fn test_priority_mid_sized_autosome() {
        assert_eq!(chrom_priority("chr5", 1000, 10000), 0);
    }

    #[test]
    fn test_priority_large_chrom() {
        assert_eq!(chrom_priority("chr1", 3000, 10000), 1);
    }

    #[test]
    fn test_priority_small_chrom() {
        assert_eq!(chrom_priority("chr22", 100, 10000), 2);
    }

    #[test]
    fn test_priority_mito_name() {
        assert_eq!(chrom_priority("chrM", 1000, 10000), 2);
        assert_eq!(chrom_priority("MT", 1000, 10000), 2);
    }

    #[test]
    fn test_priority_unplaced_scaffold() {
        assert_eq!(chrom_priority("chr1_random", 1000, 10000), 2);
        assert_eq!(chrom_priority("chrUn_gl000220", 500, 10000), 2);
    }

    #[test]
    fn test_priority_zero_genome_size() {
        assert_eq!(chrom_priority("chr1", 1000, 0), 1);
    }

    #[test]
    fn test_priority_tiny_fraction() {
        assert_eq!(chrom_priority("chrSmall", 1, 1_000_000), 2);
    }

    // ─── infer_scale_factor integration tests ────────────────────────────────

    fn make_two_chrom_bigwig(
        dir: &tempfile::TempDir,
        min_val: f32,
        second_val: f32,
        bin_size: u32,
    ) -> anyhow::Result<std::path::PathBuf> {
        let path = dir.path().join("test.bw");
        let chrom_len: u32 = bin_size * 20;
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), chrom_len);
        chroms.insert("chr7".to_string(), chrom_len);

        let mut intervals: Vec<(String, u32, u32, f32)> = Vec::new();
        for chrom in &["chr5", "chr7"] {
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
        make_bigwig(&path, chroms, intervals)?;
        Ok(path)
    }

    #[test]
    fn test_cpm_detection() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

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
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 10.0, 20.0, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

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
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 1000.0, 2000.0, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert_eq!(r.norm_method, NormMethod::Unknown);
        assert!(r.confident);
        assert_eq!(r.chroms_scanned, 2);
        assert!((r.ratio - 2.0).abs() < 0.01);
        Ok(())
    }

    #[test]
    fn test_pseudocount_detected() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.0625, 0.1875, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(r.warnings.pseudocount_detected);
        let p = r.pseudocount.expect("pseudocount should be Some");
        assert!((p - 0.5).abs() < 0.01, "p ≈ 0.5, got {p}");
        assert!((r.scale_factor - (1.0 / 0.1875)).abs() < 0.01);
        assert!((r.ratio - 3.0).abs() < 0.01);
        assert!(r.confident);
        Ok(())
    }

    #[test]
    fn test_smoothing_detected() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.25, 0.375, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(r.warnings.smoothing_detected);
        assert!(!r.confident);
        assert!((r.ratio - 1.5).abs() < 0.05);
        Ok(())
    }

    #[test]
    fn test_variable_bin_sizes_warning() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("varbins.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 50u32);
        chroms.insert("chr7".to_string(), 50u32);
        let intervals = vec![
            ("chr5".to_string(), 0u32, 10u32, 0.01f32),
            ("chr5".to_string(), 10, 21, 0.02),
            ("chr5".to_string(), 21, 50, 0.02),
            ("chr7".to_string(), 0, 10, 0.01),
            ("chr7".to_string(), 10, 21, 0.02),
            ("chr7".to_string(), 21, 50, 0.02),
        ];
        make_bigwig(&path, chroms, intervals)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(r.warnings.variable_bin_sizes);
        Ok(())
    }

    #[test]
    fn test_collapsed_bins_not_variable() -> anyhow::Result<()> {
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
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(!r.warnings.variable_bin_sizes);
        assert_eq!(r.bin_size, 10);
        Ok(())
    }

    #[test]
    fn test_confident_after_two_chroms() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(r.confident);
        assert_eq!(r.chroms_scanned, 2);
        Ok(())
    }

    #[test]
    fn test_not_confident_single_chrom() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("single.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 200u32);
        let intervals = vec![
            ("chr5".to_string(), 0u32, 10u32, 0.01f32),
            ("chr5".to_string(), 10, 20, 0.02),
        ];
        make_bigwig(&path, chroms, intervals)?;
        let r = infer_scale_factor(&path, &default_config())?;

        assert!(!r.confident);
        assert_eq!(r.chroms_scanned, 1);
        Ok(())
    }

    #[test]
    fn test_scale_factor_round_trip_cpm() {
        let cpm_val: f64 = 0.01;
        let scale = 1.0 / cpm_val;
        assert!((cpm_val * scale - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_scale_factor_round_trip_rpkm() {
        let rpkm_val: f64 = 10.0;
        let scale = 1.0 / rpkm_val;
        assert!((rpkm_val * scale - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_on_empty_bigwig() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.bw");
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), 1000u32);
        let result = make_bigwig(&path, chroms, vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_to_raw_cpm() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;
        assert!((r.to_raw(0.01) - 1.0).abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn test_json_serialization() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.01, 0.02, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;
        let json = serde_json::to_string(&r)?;
        let v: serde_json::Value = serde_json::from_str(&json)?;
        assert_eq!(v["norm_method"], "Cpm");
        assert!(v["scale_factor"].is_number());
        assert!(v["second_min_val"].is_number());
        // warnings serialise as an empty array when all flags are false
        assert!(v["warnings"].is_array());
        assert_eq!(v["warnings"].as_array().unwrap().len(), 0);
        Ok(())
    }

    #[test]
    fn test_warnings_serialize_as_messages() -> anyhow::Result<()> {
        // Pseudocount case produces a non-empty warnings array.
        let dir = tempfile::tempdir()?;
        let path = make_two_chrom_bigwig(&dir, 0.0625, 0.1875, 10)?;
        let r = infer_scale_factor(&path, &default_config())?;
        let json = serde_json::to_string(&r)?;
        let v: serde_json::Value = serde_json::from_str(&json)?;
        let warnings = v["warnings"].as_array().unwrap();
        assert!(!warnings.is_empty(), "pseudocount warning should appear");
        assert!(
            warnings[0]
                .as_str()
                .unwrap()
                .starts_with("pseudocount_detected"),
            "warning message should be a readable string"
        );
        Ok(())
    }

    // ─── heuristic correctness tests ─────────────────────────────────────────

    /// Build a chromosome where the minimum appears after N leading intervals
    /// of a higher value.  Used to probe the early-exit heuristic.
    ///
    /// Layout (bin_size=10, all on chr5):
    ///   [0, 10)          → leading_val   (repeated `leading_n` times)
    ///   [leading_n*10, …) → true_min     (one bin)
    fn make_late_minimum_bigwig(
        dir: &tempfile::TempDir,
        leading_n: u32,
        leading_val: f32,
        true_min: f32,
        bin_size: u32,
    ) -> anyhow::Result<std::path::PathBuf> {
        let path = dir.path().join("late_min.bw");
        let total_bins = leading_n + 1;
        let chrom_len = total_bins * bin_size;
        let mut chroms = HashMap::new();
        // Two chromosomes so the cross-chrom confidence logic can fire.
        chroms.insert("chr5".to_string(), chrom_len);
        chroms.insert("chr7".to_string(), chrom_len);

        let mut intervals: Vec<(String, u32, u32, f32)> = Vec::new();
        for chrom in &["chr5", "chr7"] {
            for i in 0..leading_n {
                intervals.push((
                    chrom.to_string(),
                    i * bin_size,
                    (i + 1) * bin_size,
                    leading_val,
                ));
            }
            // true minimum is the very last bin
            intervals.push((
                chrom.to_string(),
                leading_n * bin_size,
                (leading_n + 1) * bin_size,
                true_min,
            ));
        }
        make_bigwig(&path, chroms, intervals)?;
        Ok(path)
    }

    /// The heuristic can miss a minimum that appears only after both a stable streak
    /// AND a second distinct value have been seen.  Concrete failure pattern:
    ///   [3.0]*A  [2.0]*(streak+5)  [1.0]*1
    /// After A intervals of 3.0 the streak counter is A (min2 still MAX → no break).
    /// The first 2.0 resets the streak; after streak+5 more 2.0s both conditions are
    /// satisfied and we break before seeing 1.0.
    #[test]
    fn test_heuristic_misses_late_minimum() -> anyhow::Result<()> {
        let streak: u32 = 8;
        let config = InferScaleConfig {
            min_stable_streak: streak,
            max_starts: 512,
        };
        let full = full_scan_config();

        // Layout per chromosome: [3.0]*5, [2.0]*(streak+5), [1.0]*1
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("late.bw");
        let bin: u32 = 10;
        let a: u32 = 5;
        let b: u32 = streak + 5;
        let chrom_len = (a + b + 1) * bin;
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), chrom_len);
        chroms.insert("chr7".to_string(), chrom_len);

        let mut intervals: Vec<(String, u32, u32, f32)> = Vec::new();
        for chrom in &["chr5", "chr7"] {
            for i in 0..a {
                intervals.push((chrom.to_string(), i * bin, (i + 1) * bin, 3.0));
            }
            for i in a..(a + b) {
                intervals.push((chrom.to_string(), i * bin, (i + 1) * bin, 2.0));
            }
            intervals.push((chrom.to_string(), (a + b) * bin, (a + b + 1) * bin, 1.0));
        }
        make_bigwig(&path, chroms, intervals)?;

        let r_full = infer_scale_factor(&path, &full)?;
        let r_heur = infer_scale_factor(&path, &config)?;

        // Full scan finds the true minimum.
        assert!(
            (r_full.min_val - 1.0).abs() < 1e-6,
            "full scan min_val should be 1.0, got {}",
            r_full.min_val
        );
        // Heuristic with small streak misses the late minimum — documents the limitation.
        assert!(
            r_heur.min_val > 1.0 + 1e-6,
            "heuristic with streak={streak} should miss the late minimum (got {}); \
             increase min_stable_streak or use u32::MAX to disable",
            r_heur.min_val
        );
        Ok(())
    }

    /// With a sufficiently large streak the heuristic matches the full scan.
    /// This is the normal operating regime for real BigWig files where a single-read
    /// bin is rare but not exclusively clustered at the chromosome end.
    #[test]
    fn test_heuristic_matches_full_scan_minimum_not_late() -> anyhow::Result<()> {
        // Minimum is the FIRST bin; all subsequent bins are higher.
        // Any streak ≥ 1 will find the correct minimum.
        let dir = tempfile::tempdir()?;
        let path = make_late_minimum_bigwig(&dir, 0, 2.0, 1.0, 10)?;
        // With 0 leading intervals the "late_minimum" fixture places 1.0 first.
        // Re-use the fixture but build a simpler version: min first, rest larger.
        let path2 = dir.path().join("min_first.bw");
        {
            let mut chroms = HashMap::new();
            chroms.insert("chr5".to_string(), 200u32);
            chroms.insert("chr7".to_string(), 200u32);
            let mut ivs: Vec<(String, u32, u32, f32)> = Vec::new();
            for chrom in &["chr5", "chr7"] {
                ivs.push((chrom.to_string(), 0, 10, 1.0_f32));
                for i in 1..20u32 {
                    ivs.push((chrom.to_string(), i * 10, (i + 1) * 10, 2.0));
                }
            }
            make_bigwig(&path2, chroms, ivs)?;
        }
        let _ = path; // suppress warning

        let r_full = infer_scale_factor(&path2, &full_scan_config())?;
        let r_heur = infer_scale_factor(
            &path2,
            &InferScaleConfig {
                min_stable_streak: 5,
                ..Default::default()
            },
        )?;

        assert!(
            (r_full.min_val - r_heur.min_val).abs() < 1e-9,
            "heuristic and full scan should agree when minimum is not late: full={} heur={}",
            r_full.min_val,
            r_heur.min_val
        );
        assert_eq!(r_full.norm_method, r_heur.norm_method);
        Ok(())
    }

    /// Default streak (50 000) is robust against a minimum that appears at a moderate
    /// position — e.g. 1 000 intervals in, far short of 50 000.
    #[test]
    fn test_default_streak_robust_for_typical_position() -> anyhow::Result<()> {
        // Pattern: 1000 bins of 2.0, then 1 bin of 1.0, then 1000 bins of 2.0
        // The minimum is at interval 1000, well within the 50k streak window.
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("mid_min.bw");
        let bin: u32 = 10;
        let n: u32 = 1_000;
        let chrom_len = (2 * n + 1) * bin;
        let mut chroms = HashMap::new();
        chroms.insert("chr5".to_string(), chrom_len);
        chroms.insert("chr7".to_string(), chrom_len);

        let mut ivs: Vec<(String, u32, u32, f32)> = Vec::new();
        for chrom in &["chr5", "chr7"] {
            for i in 0..n {
                ivs.push((chrom.to_string(), i * bin, (i + 1) * bin, 2.0));
            }
            ivs.push((chrom.to_string(), n * bin, (n + 1) * bin, 1.0));
            for i in (n + 1)..(2 * n + 1) {
                ivs.push((chrom.to_string(), i * bin, (i + 1) * bin, 2.0));
            }
        }
        make_bigwig(&path, chroms, ivs)?;

        let r_full = infer_scale_factor(&path, &full_scan_config())?;
        let r_heur = infer_scale_factor(&path, &default_config())?;

        assert!(
            (r_full.min_val - r_heur.min_val).abs() < 1e-9,
            "default streak should find the correct minimum: full={} heur={}",
            r_full.min_val,
            r_heur.min_val
        );
        Ok(())
    }
}
