//! # Signal Normalization Module
//!
//! This module provides methods for normalizing genomic coverage signals.
//! It supports common normalization techniques used in sequencing analysis, such as:
//! *   **Raw**: No normalization.
//! *   **RPKM**: Reads Per Kilobase per Million mapped reads.
//! *   **CPM**: Counts Per Million mapped reads.

use anyhow::Result;
use log::warn;

/// Normalization methods for genomic coverage signals.
#[derive(Debug, Clone, clap::ValueEnum)]
pub enum NormalizationMethod {
    /// No normalization.
    Raw,
    /// Reads Per Kilobase per Million mapped reads.
    RPKM,
    /// Counts Per Million mapped reads.
    CPM,
}

impl NormalizationMethod {
    /// Parse a string into a `NormalizationMethod`.
    ///
    /// # Arguments
    ///
    /// * `s` - The string to parse.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Result<NormalizationMethod> {
        match s.to_lowercase().as_str() {
            "raw" => Ok(NormalizationMethod::Raw),
            "rpkm" => Ok(NormalizationMethod::RPKM),
            "cpm" => Ok(NormalizationMethod::CPM),
            _ => {
                warn!("Unknown normalization method: {s}. Defaulting to Raw");
                Ok(NormalizationMethod::Raw)
            }
        }
    }

    /// Calculate the scaling factor for the given normalization method.
    /// # Arguments:
    /// * `base_scale` - The base scaling factor (e.g., 1.0 for raw counts).
    /// * `bin_size` - The size of the genomic bin in base pairs.
    /// * `n_reads` - The total number of reads in the dataset.
    /// # Returns:
    /// A scaling factor to apply to raw counts based on the normalization method.
    /// This factor is used to adjust the raw counts to account for differences in sequencing depth and
    /// genomic feature length.
    pub fn scale_factor(&self, base_scale: f32, bin_size: u64, n_reads: u64) -> f64 {
        let base = base_scale as f64;

        match self {
            Self::Raw => base,

            Self::CPM => {
                // Counts Per Million: normalize by library size
                base * (1_000_000.0 / n_reads as f64)
            }

            Self::RPKM => {
                // Reads Per Kilobase per Million: normalize by both length and library size
                let reads_per_million = 1_000_000.0 / n_reads as f64;
                let per_kilobase = 1_000.0 / bin_size as f64; // Convert bp to kb
                base * reads_per_million * per_kilobase
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scaling_factors() {
        let bin_size = 1000; // 1kb bins
        let n_reads = 10_000_000; // 10M reads
        let base_scale = 1.0;

        let raw_factor = NormalizationMethod::Raw.scale_factor(base_scale, bin_size, n_reads);
        let cpm_factor = NormalizationMethod::CPM.scale_factor(base_scale, bin_size, n_reads);
        let rpkm_factor = NormalizationMethod::RPKM.scale_factor(base_scale, bin_size, n_reads);

        println!("Raw factor: {raw_factor}"); // Should be 1.0
        println!("CPM factor: {cpm_factor}"); // Should be 0.1 (1M/10M)
        println!("RPKM factor: {rpkm_factor}"); // Should be 0.1 (1M/10M * 1000/1K)

        // For 1000 reads in a 1kb bin with 10M total reads:
        let raw_value = 1000.0 * raw_factor; // 1000
        let cpm_value = 1000.0 * cpm_factor; // 100 CPM
        let rpkm_value = 1000.0 * rpkm_factor; // 100 RPKM

        assert_eq!(raw_value, 1000.0);
        assert_eq!(cpm_value, 100.0);
        assert_eq!(rpkm_value, 100.0); // Changed from 0.1 to 100.0
    }
}
