//! # BamNado Python Interface
//!
//! This crate provides the Python bindings for the `bamnado` library using `pyo3`.
//! It exposes high-performance Rust implementations of BAM processing tools to Python.
//!
//! ## Exposed Functions
//!
//! *   `get_signal_for_chromosome`: Calculates coverage signal for a chromosome.

use bamnado::{bam_utils, coverage_analysis, read_filter::BamReadFilter, signal_normalization};
use ndarray::prelude::*;
use numpy::{PyArray1, prelude::*};
use pyo3::exceptions::{PyKeyError, PyRuntimeError};
use pyo3::prelude::*;

fn anyhow_to_pyerr(err: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

#[pymodule]
mod _bamnado {
    use super::*;

    /// Get the signal for a specific chromosome from a BAM file.
    ///
    /// This function calculates the coverage signal for a given chromosome, optionally scaling it
    /// and using fragment length information.
    ///
    /// Args:
    ///     bam_path (str): Path to the BAM file.
    ///     chromosome_name (str): Name of the chromosome to analyze.
    ///     bin_size (int): Size of the bins for coverage calculation.
    ///     scale_factor (float): Factor to scale the signal by.
    ///     use_fragment (bool): Whether to use fragment length for coverage (True) or just read start (False).
    ///     ignore_scaffold_chromosomes (bool): Whether to ignore scaffold chromosomes during analysis.
    ///
    /// Returns:
    ///     numpy.ndarray: A 1D numpy array of floats representing the signal across the chromosome.
    ///
    /// Raises:
    ///     RuntimeError: If there is an error reading the BAM file or processing the data.
    ///     KeyError: If the specified chromosome is not found in the BAM file.
    #[pyfunction]
    #[pyo3(
        text_signature = "(bam_path, chromosome_name, bin_size, scale_factor, use_fragment, ignore_scaffold_chromosomes)"
    )]
    fn get_signal_for_chromosome(
        py: Python,
        bam_path: &str,
        chromosome_name: &str,
        bin_size: u64,
        scale_factor: f32,
        use_fragment: bool,
        ignore_scaffold_chromosomes: bool,
    ) -> PyResult<Py<PyArray1<f32>>> {
        let stats = bam_utils::BamStats::new(bam_path.into()).map_err(anyhow_to_pyerr)?;
        let chrom_size = match stats
            .chromsizes_ref_name()
            .map_err(anyhow_to_pyerr)?
            .get(chromosome_name)
        {
            Some(size) => *size,
            None => {
                return Err(PyKeyError::new_err(format!(
                    "Chromosome {} not found in BAM file.",
                    chromosome_name
                )));
            }
        };

        let bam_pileup = coverage_analysis::BamPileup::new(
            bam_path.into(),
            bin_size,
            signal_normalization::NormalizationMethod::Raw,
            scale_factor,
            use_fragment,
            BamReadFilter::default(),
            true,
            ignore_scaffold_chromosomes,
            None,
            None,
        );

        let signal = bam_pileup
            .pileup_chromosome(chromosome_name)
            .map_err(anyhow_to_pyerr)?;

        // Put the signal into an ndarray then convert to numpy array
        let mut array = Array1::zeros(chrom_size as usize);
        for interval in signal {
            let start = interval.start;
            let end = interval.stop;
            let value = interval.val as f32;
            array.slice_mut(s![start..end]).fill(value);
        }

        let numpy_array = array.into_pyarray(py).unbind();
        Ok(numpy_array)
    }

    #[pyfunction]
    fn __version__() -> &'static str {
        env!("CARGO_PKG_VERSION")
    }
}
