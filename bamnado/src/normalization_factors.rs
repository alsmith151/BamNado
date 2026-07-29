//! # Normalization Factors Module
//!
//! This module estimates between-sample scaling factors from a shared genomic bin count
//! matrix (see `bin_counts`). It supports:
//! *   **TMM** (trimmed mean of M-values), matching `edgeR::calcNormFactors(method = "TMM")`.
//! *   **csaw-background**, TMM restricted to large background bins with the most-enriched
//!     bins excluded, matching the `csaw` composition-bias workflow.
//! *   **CPM**, a naive total-count baseline for comparison.
//! *   **Median-of-ratios** (DESeq2-style), a robust alternative estimator.
//! *   **Spike-in**, using exogenous mapped-read counts read directly from the BAM index.
//!
//! ## Scale factor definition
//!
//! For the bin-count-based methods (TMM, csaw-background, median-of-ratios), each sample's
//! *effective library size* is `E_j = N_j * f_j`, where `N_j` is the sample's post-filter
//! library size and `f_j` is the method's composition-bias norm factor (geometric-mean
//! centred so that `prod(f_j) == 1`). The reported `scale_factor_j` is
//! `geometric_mean(E) / E_j`, so factors are centred on 1 and only relative depth /
//! composition differences move a sample's factor away from 1. `CPM` is a special case:
//! `scale_factor_j = 1e6 / N_j`, matching [`crate::signal_normalization::NormalizationMethod::CPM`].
//!
//! Apply the reported `scale_factors` directly to `bam-coverage -f` / `bigwig-aggregate
//! --scale-factors`.

use std::collections::HashSet;

use anyhow::{Context, Result, bail};
use ndarray::Array2;
use serde::Serialize;

use crate::bam_utils::BamStats;
use crate::bin_counts::BinCounts;

/// Between-sample normalization methods.
#[derive(Debug, Clone, clap::ValueEnum)]
pub enum NormFactorMethod {
    /// Trimmed mean of M-values (edgeR `calcNormFactors(method = "TMM")`).
    Tmm,
    /// TMM restricted to large background bins with the most-enriched bins excluded (csaw-style).
    CsawBackground,
    /// Counts per million mapped reads — a naive depth-only baseline.
    Cpm,
    /// DESeq2-style median-of-ratios estimator.
    MedianOfRatios,
    /// Spike-in (exogenous) mapped-read count ratio, read from the BAM index.
    SpikeIn,
}

/// Tuning parameters shared across the bin-count-based estimators.
#[derive(Debug, Clone)]
pub struct NormFactorParams {
    /// Sample to use as the TMM reference. `None` picks the sample whose bin total is
    /// closest to the geometric mean of all bin totals.
    pub reference_sample: Option<String>,
    /// Fraction of extreme M-values (log-ratios) trimmed from each tail before averaging.
    pub logratio_trim: f64,
    /// Fraction of extreme A-values (log-intensities) trimmed from each tail before averaging.
    pub sum_trim: f64,
    /// Minimum A-value (log-intensity) for a bin to be considered; matches edgeR's `Acutoff`.
    pub a_cutoff: f64,
    /// Weight M-values by their estimated asymptotic variance (edgeR's default behaviour).
    pub do_weighting: bool,
    /// Percentage of highest-mean-count bins dropped before estimating (csaw-background only).
    pub exclude_top_percent: f64,
}

impl Default for NormFactorParams {
    fn default() -> Self {
        Self {
            reference_sample: None,
            logratio_trim: 0.30,
            sum_trim: 0.05,
            a_cutoff: -1e10,
            do_weighting: true,
            exclude_top_percent: 5.0,
        }
    }
}

/// The result of a between-sample scaling factor computation.
#[derive(Debug, Clone, Serialize)]
pub struct ScaleFactors {
    /// Name of the method used to compute these factors.
    pub method: String,
    /// Sample names, in the same order as the other per-sample fields.
    pub sample_names: Vec<String>,
    /// Post-filter library size (or mapped-read count, for spike-in) per sample.
    pub library_sizes: Vec<u64>,
    /// Method-specific composition-bias factor per sample (geometric-mean centred to 1
    /// where applicable).
    pub norm_factors: Vec<f64>,
    /// Multiplier to pass to `bam-coverage -f` / `bigwig-aggregate --scale-factors`.
    pub scale_factors: Vec<f64>,
    /// The sample used as the TMM reference, if applicable.
    pub reference_sample: Option<String>,
    /// Total number of bins in the input count matrix.
    pub n_bins_total: usize,
    /// Number of bins actually used by the estimator (after any pre-filtering).
    pub n_bins_used: usize,
}

/// A pluggable between-sample scaling factor estimation strategy.
///
/// New methods implement this trait without changing any call sites.
pub trait NormFactorEstimator {
    /// Short, stable name for this method (used in [`ScaleFactors::method`]).
    fn name(&self) -> &'static str;

    /// Estimate scaling factors from a shared bin count matrix.
    fn estimate(&self, counts: &BinCounts, params: &NormFactorParams) -> Result<ScaleFactors>;
}

fn geometric_mean(xs: &[f64]) -> f64 {
    let logs: Vec<f64> = xs
        .iter()
        .filter(|&&x| x > 0.0 && x.is_finite())
        .map(|&x| x.ln())
        .collect();
    if logs.is_empty() {
        1.0
    } else {
        (logs.iter().sum::<f64>() / logs.len() as f64).exp()
    }
}

fn median(xs: &[f64]) -> f64 {
    let mut v: Vec<f64> = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).expect("NaN in median input"));
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Centre effective library sizes into scale factors: `geometric_mean(E) / E_j`.
fn centred_scale_factors(effective_library_sizes: &[f64]) -> Vec<f64> {
    let gm = geometric_mean(effective_library_sizes);
    effective_library_sizes
        .iter()
        .map(|&e| if e > 0.0 { gm / e } else { 1.0 })
        .collect()
}

/// Trimmed mean of M-values (edgeR `calcNormFactors(method = "TMM")`).
pub struct Tmm;

impl Tmm {
    /// Compute the pairwise TMM factor of sample `j` against reference sample `r`.
    ///
    /// Implements edgeR's rank-based double trim (ties broken by first occurrence, i.e.
    /// edgeR's `ties = "first"`), which is the one place a naive value-based trim
    /// implementation would diverge from R.
    fn pairwise_factor(
        &self,
        counts: &BinCounts,
        j: usize,
        r: usize,
        n_j: f64,
        n_r: f64,
        params: &NormFactorParams,
    ) -> f64 {
        let mut m = Vec::new();
        let mut a = Vec::new();
        let mut v = Vec::new();

        for i in 0..counts.bins.len() {
            let y_j = counts.counts[[i, j]] as f64;
            let y_r = counts.counts[[i, r]] as f64;
            if y_j <= 0.0 || y_r <= 0.0 {
                continue;
            }
            let m_g = ((y_j / n_j) / (y_r / n_r)).log2();
            let a_g = 0.5 * ((y_j / n_j) * (y_r / n_r)).log2();
            if !m_g.is_finite() || !a_g.is_finite() {
                continue;
            }
            if a_g <= params.a_cutoff {
                continue;
            }
            let v_g = (n_j - y_j) / (n_j * y_j) + (n_r - y_r) / (n_r * y_r);
            m.push(m_g);
            a.push(a_g);
            v.push(v_g);
        }

        if m.is_empty() {
            return 1.0;
        }

        let max_abs_m = m.iter().fold(0.0_f64, |acc, &x| acc.max(x.abs()));
        if max_abs_m < 1e-6 {
            return 1.0;
        }

        let n = m.len();
        let mut order_m: Vec<usize> = (0..n).collect();
        order_m.sort_by(|&x, &y| m[x].partial_cmp(&m[y]).unwrap());
        let mut rank_m = vec![0usize; n];
        for (rank, &idx) in order_m.iter().enumerate() {
            rank_m[idx] = rank;
        }

        let mut order_a: Vec<usize> = (0..n).collect();
        order_a.sort_by(|&x, &y| a[x].partial_cmp(&a[y]).unwrap());
        let mut rank_a = vec![0usize; n];
        for (rank, &idx) in order_a.iter().enumerate() {
            rank_a[idx] = rank;
        }

        // edgeR uses 1-based ranks: loL = floor(n*trim)+1; hiL = n+1-loL.
        let lo_l = ((n as f64 * params.logratio_trim).floor() as usize) + 1;
        let hi_l = n + 1 - lo_l;
        let lo_s = ((n as f64 * params.sum_trim).floor() as usize) + 1;
        let hi_s = n + 1 - lo_s;

        let mut sum_wm = 0.0;
        let mut sum_w = 0.0;
        let mut kept_m = Vec::new();

        for g in 0..n {
            let rm1 = rank_m[g] + 1;
            let ra1 = rank_a[g] + 1;
            if rm1 < lo_l || rm1 > hi_l || ra1 < lo_s || ra1 > hi_s {
                continue;
            }
            kept_m.push(m[g]);
            if params.do_weighting {
                let w = 1.0 / v[g];
                sum_wm += w * m[g];
                sum_w += w;
            }
        }

        let f = if kept_m.is_empty() {
            0.0
        } else if params.do_weighting {
            if sum_w > 0.0 { sum_wm / sum_w } else { 0.0 }
        } else {
            kept_m.iter().sum::<f64>() / kept_m.len() as f64
        };

        let f = if f.is_nan() { 0.0 } else { f };
        2f64.powf(f)
    }
}

impl NormFactorEstimator for Tmm {
    fn name(&self) -> &'static str {
        "tmm"
    }

    fn estimate(&self, counts: &BinCounts, params: &NormFactorParams) -> Result<ScaleFactors> {
        let n_samples = counts.sample_names.len();
        anyhow::ensure!(
            n_samples >= 2,
            "TMM requires at least 2 samples, got {n_samples}"
        );

        // N_j = colSums(counts) -- bin totals, matching edgeR (not post-filter library size).
        let bin_totals: Vec<f64> = (0..n_samples)
            .map(|j| {
                (0..counts.bins.len())
                    .map(|i| counts.counts[[i, j]] as f64)
                    .sum()
            })
            .collect();

        let ref_idx = match &params.reference_sample {
            Some(name) => counts
                .sample_names
                .iter()
                .position(|n| n == name)
                .ok_or_else(|| anyhow::anyhow!("Unknown TMM reference sample: {name}"))?,
            None => {
                let target = geometric_mean(&bin_totals);
                bin_totals
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        (*a - target)
                            .abs()
                            .partial_cmp(&(*b - target).abs())
                            .unwrap()
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0)
            }
        };

        let mut norm_factors = vec![1.0f64; n_samples];
        for j in 0..n_samples {
            if j == ref_idx {
                continue;
            }
            norm_factors[j] = self.pairwise_factor(
                counts,
                j,
                ref_idx,
                bin_totals[j],
                bin_totals[ref_idx],
                params,
            );
        }

        // Centre so the geometric mean of norm_factors is 1.
        let gm = geometric_mean(&norm_factors);
        for f in norm_factors.iter_mut() {
            *f /= gm;
        }

        let effective_library_sizes: Vec<f64> = counts
            .library_sizes
            .iter()
            .zip(norm_factors.iter())
            .map(|(&lib, &f)| lib as f64 * f)
            .collect();
        let scale_factors = centred_scale_factors(&effective_library_sizes);

        Ok(ScaleFactors {
            method: self.name().to_string(),
            sample_names: counts.sample_names.clone(),
            library_sizes: counts.library_sizes.clone(),
            norm_factors,
            scale_factors,
            reference_sample: Some(counts.sample_names[ref_idx].clone()),
            n_bins_total: counts.bins.len(),
            n_bins_used: counts.bins.len(),
        })
    }
}

/// TMM restricted to large background bins, with the most-enriched bins excluded.
///
/// Blacklist exclusion happens upstream (via `--blacklist` on the read filter, which zeroes
/// out those bins before they ever reach this estimator). This estimator additionally drops
/// the top `exclude_top_percent` of bins by mean count (likely enriched / peak bins) and any
/// all-zero bins, then delegates to `inner` (TMM by default).
pub struct CsawBackground {
    /// The estimator applied to the filtered background bins.
    pub inner: Box<dyn NormFactorEstimator>,
}

impl Default for CsawBackground {
    fn default() -> Self {
        Self {
            inner: Box::new(Tmm),
        }
    }
}

impl NormFactorEstimator for CsawBackground {
    fn name(&self) -> &'static str {
        "csaw-background"
    }

    fn estimate(&self, counts: &BinCounts, params: &NormFactorParams) -> Result<ScaleFactors> {
        let n_bins_total = counts.bins.len();
        let n_samples = counts.sample_names.len();

        let mut mean_counts: Vec<(usize, f64)> = (0..n_bins_total)
            .map(|i| {
                let mean = (0..n_samples)
                    .map(|j| counts.counts[[i, j]] as f64)
                    .sum::<f64>()
                    / n_samples as f64;
                (i, mean)
            })
            .collect();

        let n_exclude =
            ((n_bins_total as f64) * (params.exclude_top_percent / 100.0)).round() as usize;
        mean_counts.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        let excluded: HashSet<usize> = mean_counts
            .iter()
            .take(n_exclude.min(n_bins_total))
            .map(|(i, _)| *i)
            .collect();

        let mut keep_rows = Vec::new();
        let mut keep_bins = Vec::new();
        for i in 0..n_bins_total {
            if excluded.contains(&i) {
                continue;
            }
            let all_zero = (0..n_samples).all(|j| counts.counts[[i, j]] == 0);
            if all_zero {
                continue;
            }
            keep_rows.push(i);
            keep_bins.push(counts.bins[i].clone());
        }

        anyhow::ensure!(
            !keep_rows.is_empty(),
            "csaw-background: no bins remain after excluding the top {}% and all-zero bins",
            params.exclude_top_percent
        );

        let mut filtered_counts = Array2::<u64>::zeros((keep_rows.len(), n_samples));
        for (new_i, &old_i) in keep_rows.iter().enumerate() {
            for j in 0..n_samples {
                filtered_counts[[new_i, j]] = counts.counts[[old_i, j]];
            }
        }

        let filtered = BinCounts {
            sample_names: counts.sample_names.clone(),
            bin_size: counts.bin_size,
            bins: keep_bins,
            counts: filtered_counts,
            library_sizes: counts.library_sizes.clone(),
        };

        let mut result = self
            .inner
            .estimate(&filtered, params)
            .context("csaw-background: inner estimator failed")?;
        result.method = self.name().to_string();
        result.n_bins_total = n_bins_total;
        result.n_bins_used = keep_rows.len();
        Ok(result)
    }
}

/// Counts per million mapped reads — a naive depth-only baseline for comparison.
///
/// Matches [`crate::signal_normalization::NormalizationMethod::CPM`]: `scale_factor_j = 1e6 /
/// library_size_j`. `norm_factors` are trivially 1.0 (no composition-bias correction).
pub struct Cpm;

impl NormFactorEstimator for Cpm {
    fn name(&self) -> &'static str {
        "cpm"
    }

    fn estimate(&self, counts: &BinCounts, _params: &NormFactorParams) -> Result<ScaleFactors> {
        let n_samples = counts.sample_names.len();
        let norm_factors = vec![1.0; n_samples];
        let scale_factors: Vec<f64> = counts
            .library_sizes
            .iter()
            .map(|&lib| if lib == 0 { 0.0 } else { 1e6 / lib as f64 })
            .collect();

        Ok(ScaleFactors {
            method: self.name().to_string(),
            sample_names: counts.sample_names.clone(),
            library_sizes: counts.library_sizes.clone(),
            norm_factors,
            scale_factors,
            reference_sample: None,
            n_bins_total: counts.bins.len(),
            n_bins_used: counts.bins.len(),
        })
    }
}

/// DESeq2-style median-of-ratios estimator.
///
/// Uses bins where all samples have a nonzero count: `gm_g = exp(mean_j(ln y_gj))`,
/// `s_j = median_g(y_gj / gm_g)`, `scale_factor_j = median(s) / s_j`.
pub struct MedianOfRatios;

impl NormFactorEstimator for MedianOfRatios {
    fn name(&self) -> &'static str {
        "median-of-ratios"
    }

    fn estimate(&self, counts: &BinCounts, _params: &NormFactorParams) -> Result<ScaleFactors> {
        let n_samples = counts.sample_names.len();
        let mut per_sample_ratios: Vec<Vec<f64>> = vec![Vec::new(); n_samples];

        for i in 0..counts.bins.len() {
            let row: Vec<u64> = (0..n_samples).map(|j| counts.counts[[i, j]]).collect();
            if row.contains(&0) {
                continue;
            }
            let log_gm = row.iter().map(|&c| (c as f64).ln()).sum::<f64>() / n_samples as f64;
            let gm = log_gm.exp();
            if gm <= 0.0 || !gm.is_finite() {
                continue;
            }
            for (j, ratios) in per_sample_ratios.iter_mut().enumerate() {
                ratios.push(row[j] as f64 / gm);
            }
        }

        anyhow::ensure!(
            !per_sample_ratios[0].is_empty(),
            "median-of-ratios: no bins with a nonzero count in every sample"
        );

        let s: Vec<f64> = per_sample_ratios
            .iter()
            .map(|ratios| median(ratios))
            .collect();
        let median_s = median(&s);

        let scale_factors: Vec<f64> = s
            .iter()
            .map(|&sj| if sj > 0.0 { median_s / sj } else { 1.0 })
            .collect();
        let norm_factors: Vec<f64> = s
            .iter()
            .map(|&sj| if median_s > 0.0 { sj / median_s } else { 1.0 })
            .collect();

        Ok(ScaleFactors {
            method: self.name().to_string(),
            sample_names: counts.sample_names.clone(),
            library_sizes: counts.library_sizes.clone(),
            norm_factors,
            scale_factors,
            reference_sample: None,
            n_bins_total: counts.bins.len(),
            n_bins_used: per_sample_ratios[0].len(),
        })
    }
}

/// Return the estimator implementation for a bin-count-based [`NormFactorMethod`].
///
/// # Errors
///
/// Returns an error for [`NormFactorMethod::SpikeIn`], which does not use the bin count
/// matrix at all — use [`spike_in_scale_factors`] instead.
pub fn estimator_for(method: &NormFactorMethod) -> Result<Box<dyn NormFactorEstimator>> {
    Ok(match method {
        NormFactorMethod::Tmm => Box::new(Tmm),
        NormFactorMethod::CsawBackground => Box::new(CsawBackground::default()),
        NormFactorMethod::Cpm => Box::new(Cpm),
        NormFactorMethod::MedianOfRatios => Box::new(MedianOfRatios),
        NormFactorMethod::SpikeIn => bail!(
            "SpikeIn normalization does not use the bin count matrix; call \
             `spike_in_scale_factors` with per-sample BamStats and an exogenous-name prefix \
             instead."
        ),
    })
}

/// Compute between-sample scaling factors from a shared bin count matrix.
///
/// # Arguments
///
/// * `counts` - The shared sample × bin count matrix (see `bin_counts::count_bins`).
/// * `method` - The normalization method to use. Must not be [`NormFactorMethod::SpikeIn`]
///   (use [`spike_in_scale_factors`] for that method instead).
/// * `params` - Tuning parameters for the bin-count-based estimators.
pub fn compute_scale_factors(
    counts: &BinCounts,
    method: &NormFactorMethod,
    params: &NormFactorParams,
) -> Result<ScaleFactors> {
    estimator_for(method)?.estimate(counts, params)
}

/// Compute between-sample scaling factors from spike-in (exogenous) mapped-read counts.
///
/// Counts come directly from each BAM's index (`BamStats::mapped_reads_by_prefix`) — no BAM
/// record scan, and no read filtering (MAPQ, duplicates, etc. is not applied). This is the
/// same basis as `BamStats::n_mapped` and the existing `SplitStats::spikein_norm_factor`.
///
/// # Arguments
///
/// * `bam_stats` - Per-sample `BamStats`, same length and order as `sample_names`.
/// * `sample_names` - Sample names.
/// * `exogenous_prefix` - Reference-name prefix identifying spike-in sequences.
pub fn spike_in_scale_factors(
    bam_stats: &[BamStats],
    sample_names: &[String],
    exogenous_prefix: &str,
) -> Result<ScaleFactors> {
    anyhow::ensure!(
        !bam_stats.is_empty(),
        "At least one BAM file is required for spike-in normalization"
    );
    anyhow::ensure!(
        bam_stats.len() == sample_names.len(),
        "Number of BAM stats ({}) does not match number of sample names ({})",
        bam_stats.len(),
        sample_names.len()
    );

    let exogenous_counts: Vec<u64> = bam_stats
        .iter()
        .map(|stats| stats.mapped_reads_by_prefix(exogenous_prefix).0)
        .collect();

    for (name, &count) in sample_names.iter().zip(exogenous_counts.iter()) {
        anyhow::ensure!(
            count > 0,
            "Sample '{name}' has zero reads mapped to references with prefix '{exogenous_prefix}'"
        );
    }

    let exogenous_f: Vec<f64> = exogenous_counts.iter().map(|&c| c as f64).collect();
    let gm = geometric_mean(&exogenous_f);
    let scale_factors: Vec<f64> = exogenous_f.iter().map(|&c| gm / c).collect();
    let norm_factors: Vec<f64> = exogenous_f.iter().map(|&c| c / gm).collect();
    let library_sizes: Vec<u64> = bam_stats.iter().map(|stats| stats.n_mapped()).collect();

    Ok(ScaleFactors {
        method: "spike-in".to_string(),
        sample_names: sample_names.to_vec(),
        library_sizes,
        norm_factors,
        scale_factors,
        reference_sample: None,
        n_bins_total: 0,
        n_bins_used: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bin_counts::BinRegion;

    fn make_bins(n: usize) -> Vec<BinRegion> {
        (0..n)
            .map(|i| BinRegion {
                chrom: "chr1".to_string(),
                start: (i * 1000 + 1) as u64,
                end: ((i + 1) * 1000) as u64,
            })
            .collect()
    }

    fn make_counts(
        sample_names: Vec<&str>,
        library_sizes: Vec<u64>,
        rows: Vec<Vec<u64>>,
    ) -> BinCounts {
        let n_bins = rows.len();
        let n_samples = sample_names.len();
        let mut counts = Array2::<u64>::zeros((n_bins, n_samples));
        for (i, row) in rows.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                counts[[i, j]] = v;
            }
        }
        BinCounts {
            sample_names: sample_names.into_iter().map(String::from).collect(),
            bin_size: 1000,
            bins: make_bins(n_bins),
            counts,
            library_sizes,
        }
    }

    #[test]
    fn test_tmm_identical_samples_gives_factor_one() {
        let rows: Vec<Vec<u64>> = (0..20).map(|i| vec![10 + i, 10 + i, 10 + i]).collect();
        let counts = make_counts(vec!["a", "b", "c"], vec![1000, 1000, 1000], rows);

        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::Tmm,
            &NormFactorParams::default(),
        )
        .expect("TMM failed");

        for &f in &result.norm_factors {
            assert!((f - 1.0).abs() < 1e-9, "norm_factor {f} not ~1.0");
        }
        for &f in &result.scale_factors {
            assert!((f - 1.0).abs() < 1e-9, "scale_factor {f} not ~1.0");
        }
    }

    #[test]
    fn test_tmm_prod_norm_factors_is_one() {
        // Sample b has a handful of bins massively inflated (composition bias); a and c
        // are otherwise identical background samples.
        let mut rows: Vec<Vec<u64>> = (0..20).map(|i| vec![50 + i, 50 + i, 50 + i]).collect();
        for row in rows.iter_mut().take(3) {
            row[1] *= 50;
        }
        let counts = make_counts(vec!["a", "b", "c"], vec![2000, 3500, 2000], rows);

        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::Tmm,
            &NormFactorParams::default(),
        )
        .expect("TMM failed");

        let prod: f64 = result.norm_factors.iter().product();
        assert!((prod - 1.0).abs() < 1e-6, "prod(norm_factors) = {prod}");

        // The unbiased pair (a, c) should end up close to factor 1 relative to each other.
        let ratio = result.norm_factors[0] / result.norm_factors[2];
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "unbiased pair (a, c) diverged: ratio = {ratio}"
        );
    }

    #[test]
    fn test_tmm_reference_sample_override() {
        let rows: Vec<Vec<u64>> = (0..20).map(|i| vec![10 + i, 20 + i, 30 + i]).collect();
        let counts = make_counts(vec!["a", "b", "c"], vec![1000, 2000, 3000], rows);

        let params = NormFactorParams {
            reference_sample: Some("c".to_string()),
            ..Default::default()
        };
        let result =
            compute_scale_factors(&counts, &NormFactorMethod::Tmm, &params).expect("TMM failed");
        assert_eq!(result.reference_sample.as_deref(), Some("c"));
    }

    #[test]
    fn test_tmm_unknown_reference_sample_errors() {
        let rows: Vec<Vec<u64>> = (0..20).map(|i| vec![10 + i, 20 + i]).collect();
        let counts = make_counts(vec!["a", "b"], vec![1000, 2000], rows);

        let params = NormFactorParams {
            reference_sample: Some("nonexistent".to_string()),
            ..Default::default()
        };
        let result = compute_scale_factors(&counts, &NormFactorMethod::Tmm, &params);
        assert!(result.is_err());
    }

    #[test]
    fn test_tmm_single_sample_errors() {
        let rows: Vec<Vec<u64>> = (0..5).map(|i| vec![10 + i]).collect();
        let counts = make_counts(vec!["a"], vec![1000], rows);
        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::Tmm,
            &NormFactorParams::default(),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_cpm_trivial() {
        let rows: Vec<Vec<u64>> = (0..5).map(|i| vec![10 + i, 20 + i]).collect();
        let counts = make_counts(vec!["a", "b"], vec![1_000_000, 2_000_000], rows);
        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::Cpm,
            &NormFactorParams::default(),
        )
        .expect("CPM failed");
        assert_eq!(result.norm_factors, vec![1.0, 1.0]);
        assert!((result.scale_factors[0] - 1.0).abs() < 1e-9);
        assert!((result.scale_factors[1] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_median_of_ratios_known_depth_difference() {
        // Sample b has exactly 2x the counts of sample a in every bin.
        let rows: Vec<Vec<u64>> = (0..20).map(|i| vec![100 + i, 2 * (100 + i)]).collect();
        let counts = make_counts(vec!["a", "b"], vec![10_000, 20_000], rows);
        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::MedianOfRatios,
            &NormFactorParams::default(),
        )
        .expect("median-of-ratios failed");

        // b is twice as deep as a, so b's scale factor should be half of a's.
        let ratio = result.scale_factors[1] / result.scale_factors[0];
        assert!(
            (ratio - 0.5).abs() < 1e-6,
            "expected b's scale factor to be half of a's, ratio = {ratio}"
        );
    }

    #[test]
    fn test_csaw_background_excludes_top_bins_and_all_zero_bins() {
        let mut rows: Vec<Vec<u64>> = (0..20).map(|_| vec![10, 10]).collect();
        // Two "peak" bins with a huge count, and a couple of all-zero bins.
        rows[0] = vec![10_000, 10_000];
        rows[1] = vec![9_000, 9_000];
        rows[2] = vec![0, 0];
        rows[3] = vec![0, 0];
        let counts = make_counts(vec!["a", "b"], vec![1000, 1000], rows);

        let params = NormFactorParams {
            exclude_top_percent: 10.0,
            ..Default::default()
        };
        let result = compute_scale_factors(&counts, &NormFactorMethod::CsawBackground, &params)
            .expect("csaw-background failed");

        assert_eq!(result.n_bins_total, 20);
        // 10% of 20 = 2 excluded as top bins, plus 2 all-zero bins dropped.
        assert_eq!(result.n_bins_used, 16);
    }

    #[test]
    fn test_all_zero_bins_dropped_from_median_of_ratios() {
        let mut rows: Vec<Vec<u64>> = (0..10).map(|i| vec![10 + i, 10 + i]).collect();
        rows.push(vec![0, 0]);
        let counts = make_counts(vec!["a", "b"], vec![1000, 1000], rows);
        let result = compute_scale_factors(
            &counts,
            &NormFactorMethod::MedianOfRatios,
            &NormFactorParams::default(),
        )
        .expect("median-of-ratios failed");
        assert_eq!(result.n_bins_used, 10);
    }

    #[test]
    fn test_estimator_for_spike_in_errors() {
        let result = estimator_for(&NormFactorMethod::SpikeIn);
        assert!(result.is_err());
    }
}
