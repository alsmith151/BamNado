//! # BamNado Python Interface
//!
//! This crate provides the Python bindings for the `bamnado` library using `pyo3`.
//! It exposes high-performance Rust implementations of BAM processing tools to Python.
//!
//! ## Exposed Classes / Functions
//!
//! *   `ReadFilter`: Configuration for filtering BAM reads.
//! *   `get_signal_for_chromosome`: Calculates coverage signal for a chromosome.

use ahash::HashSet;
use bamnado::{bam_utils, coverage_analysis, read_filter::BamReadFilter, signal_normalization};
use bio_types::strand::Strand;
use ndarray::prelude::*;
use numpy::{PyArray1, prelude::*};
use pyo3::exceptions::{PyKeyError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

fn anyhow_to_pyerr(err: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

#[pymodule]
mod _bamnado {
    use super::*;

    /// Configuration for filtering BAM reads.
    ///
    /// Controls which reads are included in coverage calculations based on mapping
    /// quality, read length, strand, paired-end flags, barcodes, and genomic blacklists.
    ///
    /// Args:
    ///     min_mapq (int): Minimum mapping quality score. Default: 0.
    ///     proper_pair (bool): Keep only properly paired reads. Default: True.
    ///     min_length (int): Minimum read length in bp. Default: 0.
    ///     max_length (int): Maximum read length in bp. Default: 1000.
    ///     strand (str): Strand to keep: ``"forward"`` (``"fwd"``, ``"+"``), ``"reverse"`` (``"rev"``, ``"-"``), or ``"both"``. Default: ``"both"``.
    ///     min_fragment_length (int | None): Minimum fragment (insert) length for paired-end reads. Default: None.
    ///     max_fragment_length (int | None): Maximum fragment (insert) length for paired-end reads. Default: None.
    ///     blacklist_bed (str | None): Path to a BED file of regions to exclude. Default: None.
    ///     whitelisted_barcodes (list[str] | None): Cell barcodes (CB tag) to include. Default: None.
    ///     read_group (str | None): Read group (RG tag) to keep. Default: None.
    ///     filter_tag (str | None): Two-character SAM tag name to filter on (e.g. ``"CB"``). Default: None.
    ///     filter_tag_value (str | None): Required string value for *filter_tag*. Default: None.
    #[pyclass(from_py_object)]
    #[derive(Clone)]
    pub struct ReadFilter {
        pub min_mapq: u8,
        pub proper_pair: bool,
        pub min_length: u32,
        pub max_length: u32,
        pub strand: String,
        pub min_fragment_length: Option<u32>,
        pub max_fragment_length: Option<u32>,
        pub blacklist_bed: Option<String>,
        pub whitelisted_barcodes: Option<Vec<String>>,
        pub read_group: Option<String>,
        pub filter_tag: Option<String>,
        pub filter_tag_value: Option<String>,
    }

    #[pymethods]
    impl ReadFilter {
        #[new]
        #[pyo3(signature = (
            min_mapq=0,
            proper_pair=true,
            min_length=0,
            max_length=1000,
            strand="both",
            min_fragment_length=None,
            max_fragment_length=None,
            blacklist_bed=None,
            whitelisted_barcodes=None,
            read_group=None,
            filter_tag=None,
            filter_tag_value=None
        ))]
        #[allow(clippy::too_many_arguments)]
        fn new(
            min_mapq: u8,
            proper_pair: bool,
            min_length: u32,
            max_length: u32,
            strand: &str,
            min_fragment_length: Option<u32>,
            max_fragment_length: Option<u32>,
            blacklist_bed: Option<&str>,
            whitelisted_barcodes: Option<Vec<String>>,
            read_group: Option<String>,
            filter_tag: Option<String>,
            filter_tag_value: Option<String>,
        ) -> PyResult<Self> {
            match strand {
                "forward" | "fwd" | "+" | "reverse" | "rev" | "-" | "both" => {}
                other => {
                    return Err(PyValueError::new_err(format!(
                        "Invalid strand '{other}': must be \"forward\", \"reverse\", or \"both\""
                    )));
                }
            }
            Ok(Self {
                min_mapq,
                proper_pair,
                min_length,
                max_length,
                strand: strand.to_string(),
                min_fragment_length,
                max_fragment_length,
                blacklist_bed: blacklist_bed.map(String::from),
                whitelisted_barcodes,
                read_group,
                filter_tag,
                filter_tag_value,
            })
        }

        fn __repr__(&self) -> String {
            format!(
                "ReadFilter(min_mapq={}, proper_pair={}, min_length={}, max_length={}, strand={:?}, \
                 min_fragment_length={:?}, max_fragment_length={:?}, blacklist_bed={:?}, \
                 whitelisted_barcodes={}, read_group={:?}, filter_tag={:?}, filter_tag_value={:?})",
                self.min_mapq,
                self.proper_pair,
                self.min_length,
                self.max_length,
                self.strand,
                self.min_fragment_length,
                self.max_fragment_length,
                self.blacklist_bed,
                self.whitelisted_barcodes
                    .as_ref()
                    .map(|v| format!("<{} barcodes>", v.len()))
                    .unwrap_or_else(|| "None".to_string()),
                self.read_group,
                self.filter_tag,
                self.filter_tag_value,
            )
        }
    }

    /// Convert a `ReadFilter` pyclass into a `BamReadFilter`, loading any BED blacklist
    /// and resolving chromosome names to IDs using the provided `BamStats`.
    fn build_bam_read_filter(
        read_filter: Option<ReadFilter>,
        bam_stats: &bam_utils::BamStats,
    ) -> PyResult<BamReadFilter> {
        let rf = match read_filter {
            None => return Ok(BamReadFilter::default()),
            Some(rf) => rf,
        };

        let strand = match rf.strand.as_str() {
            "forward" | "fwd" | "+" => Strand::Forward,
            "reverse" | "rev" | "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };

        if (rf.min_fragment_length.is_some() || rf.max_fragment_length.is_some())
            && !bam_stats.is_paired_end().map_err(anyhow_to_pyerr)?
        {
            return Err(PyValueError::new_err(
                "Fragment length filtering requires paired-end reads, but the BAM file \
                 does not appear to contain paired-end data.",
            ));
        }

        let blacklisted_locations = if let Some(bed_path) = &rf.blacklist_bed {
            let lapper =
                bam_utils::bed_to_lapper(PathBuf::from(bed_path)).map_err(anyhow_to_pyerr)?;
            let lapper_ids = bam_utils::convert_lapper_chrom_names_to_ids(lapper, bam_stats)
                .map_err(anyhow_to_pyerr)?;
            Some(lapper_ids)
        } else {
            None
        };

        let whitelisted_barcodes = rf
            .whitelisted_barcodes
            .map(|barcodes| barcodes.into_iter().collect::<HashSet<String>>());

        Ok(BamReadFilter::new(
            strand,
            rf.proper_pair,
            Some(rf.min_mapq),
            Some(rf.min_length),
            Some(rf.max_length),
            rf.read_group,
            blacklisted_locations,
            whitelisted_barcodes,
            rf.filter_tag,
            rf.filter_tag_value,
            rf.min_fragment_length,
            rf.max_fragment_length,
        ))
    }

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
    ///     read_filter (ReadFilter | None): Read filter configuration. Uses default settings if None.
    ///
    /// Returns:
    ///     numpy.ndarray: A 1D numpy array of floats representing the signal across the chromosome.
    ///
    /// Raises:
    ///     RuntimeError: If there is an error reading the BAM file or processing the data.
    ///     KeyError: If the specified chromosome is not found in the BAM file.
    ///     ValueError: If an invalid strand or fragment length configuration is provided.
    #[allow(clippy::too_many_arguments)]
    #[pyfunction]
    #[pyo3(
        text_signature = "(bam_path, chromosome_name, bin_size, scale_factor, use_fragment, ignore_scaffold_chromosomes, read_filter=None)",
        signature = (bam_path, chromosome_name, bin_size, scale_factor, use_fragment, ignore_scaffold_chromosomes, read_filter=None)
    )]
    fn get_signal_for_chromosome(
        py: Python,
        bam_path: &str,
        chromosome_name: &str,
        bin_size: u64,
        scale_factor: f32,
        use_fragment: bool,
        ignore_scaffold_chromosomes: bool,
        read_filter: Option<ReadFilter>,
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

        let filter = build_bam_read_filter(read_filter, &stats)?;

        let bam_pileup = coverage_analysis::BamPileup::new(
            bam_path.into(),
            bin_size,
            signal_normalization::NormalizationMethod::Raw,
            scale_factor,
            use_fragment,
            filter,
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
