use anyhow;
use ndarray::prelude::*;
use numpy::{PyArray1, prelude::*};
use pyo3::exceptions::{PyKeyError, PyRuntimeError};
use pyo3::prelude::*;
use bamnado::{bam_utils, coverage_analysis, signal_normalization, read_filter::BamReadFilter};

fn anyhow_to_pyerr(err: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

#[pymodule]
mod _bamnado {
    use super::*;

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
        let stats = bam_utils::BamStats::new(bam_path.into()).map_err(|e| anyhow_to_pyerr(e))?;
        let chrom_size = match stats
            .chromsizes_ref_name()
            .map_err(|x| anyhow_to_pyerr(x))?
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
            .map_err(|e| anyhow_to_pyerr(e))?;

        // Put the signal into an ndarray then convert to numpy array
        let mut array = Array1::zeros(chrom_size as usize);
        for interval in signal {
            let start = interval.start as usize;
            let end = interval.stop as usize;
            let value = interval.val as f32;
            array.slice_mut(s![start..end]).fill(value);
        }

        let numpy_array = array.into_pyarray(py).unbind();
        Ok(numpy_array)
    }
}
