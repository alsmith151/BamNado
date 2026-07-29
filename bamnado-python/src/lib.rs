//! # BamNado Python Interface
//!
//! This crate provides the Python bindings for the `bamnado` library using `pyo3`.
//! It exposes high-performance Rust implementations of BAM processing tools to Python.
//!
//! ## Exposed Classes / Functions
//!
//! *   `ReadFilter`: Configuration for filtering BAM reads.
//! *   `get_signal_for_chromosome`: Calculates coverage signal for a chromosome.
//! *   `compute_scale_factors`: Computes between-sample scaling factors from BAM files.
//! *   `ScaleFactorResult`: Result of a between-sample scaling factor computation.

use ahash::HashSet;
use bamnado::{bam_utils, coverage_analysis, read_filter::BamReadFilter, signal_normalization};
use bio_types::strand::Strand;
use ndarray::prelude::*;
use numpy::{PyArray1, PyArray2, prelude::*};
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
    #[derive(Debug, Clone)]
    pub struct ReadFilter {
        #[pyo3(get, set)]
        pub min_mapq: u8,
        #[pyo3(get, set)]
        pub proper_pair: bool,
        #[pyo3(get, set)]
        pub min_length: u32,
        #[pyo3(get, set)]
        pub max_length: u32,
        #[pyo3(get, set)]
        pub strand: String,
        #[pyo3(get, set)]
        pub min_fragment_length: Option<u32>,
        #[pyo3(get, set)]
        pub max_fragment_length: Option<u32>,
        #[pyo3(get, set)]
        pub blacklist_bed: Option<String>,
        #[pyo3(get, set)]
        pub whitelisted_barcodes: Option<Vec<String>>,
        #[pyo3(get, set)]
        pub read_group: Option<String>,
        #[pyo3(get, set)]
        pub filter_tag: Option<String>,
        #[pyo3(get, set)]
        pub filter_tag_value: Option<String>,
    }

    #[pymethods]
    impl ReadFilter {
        #[new]
        #[pyo3(signature = (
            min_mapq=0,
            proper_pair=false,
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

        fn __str__(&self) -> String {
            self.__repr__()
        }

        fn __copy__(&self) -> Self {
            self.clone()
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
            let value = interval.val as f32 * scale_factor;
            array.slice_mut(s![start..end]).fill(value);
        }

        let numpy_array = array.into_pyarray(py).unbind();
        Ok(numpy_array)
    }

    /// Result of a between-sample scaling factor computation.
    ///
    /// See `compute_scale_factors` for how these fields are derived.
    #[pyclass]
    pub struct ScaleFactorResult {
        #[pyo3(get)]
        pub method: String,
        #[pyo3(get)]
        pub sample_names: Vec<String>,
        #[pyo3(get)]
        pub library_sizes: Vec<u64>,
        #[pyo3(get)]
        pub norm_factors: Vec<f64>,
        #[pyo3(get)]
        pub scale_factors: Vec<f64>,
        #[pyo3(get)]
        pub reference_sample: Option<String>,
        #[pyo3(get)]
        pub n_bins_total: usize,
        #[pyo3(get)]
        pub n_bins_used: usize,
        counts: Option<Array2<u64>>,
        bin_chroms: Vec<String>,
        bin_starts: Vec<u64>,
        bin_ends: Vec<u64>,
    }

    #[pymethods]
    impl ScaleFactorResult {
        /// The `(n_bins, n_samples)` bin count matrix, or `None` if `return_counts=False`.
        fn counts<'py>(&self, py: Python<'py>) -> Option<Py<PyArray2<u64>>> {
            self.counts.as_ref().map(|c| c.to_pyarray(py).unbind())
        }

        /// Chromosome name for each row of `counts()`. Empty if `return_counts=False`.
        #[getter]
        fn bin_chroms(&self) -> Vec<String> {
            self.bin_chroms.clone()
        }

        /// 1-based inclusive start position for each row of `counts()`.
        #[getter]
        fn bin_starts(&self) -> Vec<u64> {
            self.bin_starts.clone()
        }

        /// 1-based inclusive end position for each row of `counts()`.
        #[getter]
        fn bin_ends(&self) -> Vec<u64> {
            self.bin_ends.clone()
        }

        fn __repr__(&self) -> String {
            format!(
                "ScaleFactorResult(method={:?}, sample_names={:?}, scale_factors={:?})",
                self.method, self.sample_names, self.scale_factors
            )
        }
    }

    /// Compute between-sample scaling factors from BAM files.
    ///
    /// Bins the genome, counts filtered reads (or fragments) per bin per sample, and
    /// estimates a scaling factor for each sample that removes composition and technical
    /// bias — e.g. so that a stronger pulldown, which consumes more of its library on peaks
    /// and deflates its background, does not make every other region in that sample look
    /// artificially low relative to the others. Apply the returned `scale_factors` to
    /// ``bam-coverage -f`` / ``bigwig-aggregate --scale-factors``.
    ///
    /// Args:
    ///     bam_paths (list[str]): Input BAM file paths, one per sample.
    ///     method (str): Normalization method: ``"tmm"``, ``"csaw-background"`` (default),
    ///         ``"cpm"``, ``"median-of-ratios"``, or ``"spike-in"``.
    ///     bin_size (int): Genomic bin width in base pairs. Ignored for ``"spike-in"``.
    ///     sample_names (list[str] | None): Sample names, same length/order as *bam_paths*.
    ///         Defaults to each BAM file's stem, disambiguated if duplicated.
    ///     exclude_top_percent (float): Percentage of highest-mean-count bins dropped before
    ///         estimating (``"csaw-background"`` only). Default: 5.0.
    ///     reference_sample (str | None): Sample name to use as the TMM reference. Defaults
    ///         to the sample whose bin total is closest to the geometric mean.
    ///     logratio_trim (float): TMM trim fraction for M-values. Default: 0.3.
    ///     sum_trim (float): TMM trim fraction for A-values. Default: 0.05.
    ///     exogenous_prefix (str | None): Reference-name prefix identifying spike-in
    ///         sequences. Required for ``method="spike-in"``.
    ///     use_fragment (bool): Count fragments (pairs) instead of individual read
    ///         alignments. Default: False.
    ///     ignore_scaffold_chromosomes (bool): Skip scaffold / unplaced chromosomes when
    ///         building the bin grid. Default: True.
    ///     read_filter (ReadFilter | None): Read filter configuration, applied per sample.
    ///         Uses default settings if None.
    ///     return_counts (bool): If True, populate ``ScaleFactorResult.counts()`` with the
    ///         underlying bin count matrix (for QC / M-A plots). Default: False.
    ///
    /// Returns:
    ///     ScaleFactorResult: The computed scaling factors and (optionally) count matrix.
    ///
    /// Raises:
    ///     ValueError: If *method* or *sample_names* is invalid, or *exogenous_prefix* is
    ///         missing for ``method="spike-in"``.
    ///     RuntimeError: If a BAM file cannot be read or the count/estimation step fails
    ///         (e.g. mismatched chromosome sets between samples).
    #[allow(clippy::too_many_arguments)]
    #[pyfunction]
    #[pyo3(signature = (
        bam_paths,
        method="csaw-background",
        bin_size=10000,
        sample_names=None,
        exclude_top_percent=5.0,
        reference_sample=None,
        logratio_trim=0.3,
        sum_trim=0.05,
        exogenous_prefix=None,
        use_fragment=false,
        ignore_scaffold_chromosomes=true,
        read_filter=None,
        return_counts=false
    ))]
    fn compute_scale_factors(
        bam_paths: Vec<String>,
        method: &str,
        bin_size: u64,
        sample_names: Option<Vec<String>>,
        exclude_top_percent: f64,
        reference_sample: Option<String>,
        logratio_trim: f64,
        sum_trim: f64,
        exogenous_prefix: Option<String>,
        use_fragment: bool,
        ignore_scaffold_chromosomes: bool,
        read_filter: Option<ReadFilter>,
        return_counts: bool,
    ) -> PyResult<ScaleFactorResult> {
        use bamnado::normalization_factors::NormFactorMethod;

        let method_enum = match method {
            "tmm" => NormFactorMethod::Tmm,
            "csaw-background" => NormFactorMethod::CsawBackground,
            "cpm" => NormFactorMethod::Cpm,
            "median-of-ratios" => NormFactorMethod::MedianOfRatios,
            "spike-in" => NormFactorMethod::SpikeIn,
            other => {
                return Err(PyValueError::new_err(format!(
                    "Invalid method '{other}': must be one of \"tmm\", \"csaw-background\", \
                     \"cpm\", \"median-of-ratios\", \"spike-in\""
                )));
            }
        };

        if bam_paths.is_empty() {
            return Err(PyValueError::new_err("At least one BAM file is required"));
        }
        let bam_paths: Vec<PathBuf> = bam_paths.into_iter().map(PathBuf::from).collect();

        let sample_names: Vec<String> = match sample_names {
            Some(names) => {
                if names.len() != bam_paths.len() {
                    return Err(PyValueError::new_err(format!(
                        "sample_names length ({}) does not match bam_paths length ({})",
                        names.len(),
                        bam_paths.len()
                    )));
                }
                names
            }
            None => bam_paths
                .iter()
                .map(|p| {
                    p.file_stem()
                        .and_then(|s| s.to_str())
                        .map(String::from)
                        .unwrap_or_else(|| p.display().to_string())
                })
                .collect(),
        };

        let bam_stats: Vec<bam_utils::BamStats> = bam_paths
            .iter()
            .map(|p| bam_utils::BamStats::new(p.clone()).map_err(anyhow_to_pyerr))
            .collect::<PyResult<Vec<_>>>()?;

        let mut counts_data: Option<(Array2<u64>, Vec<String>, Vec<u64>, Vec<u64>)> = None;

        let result = if matches!(method_enum, NormFactorMethod::SpikeIn) {
            let prefix = exogenous_prefix.ok_or_else(|| {
                PyValueError::new_err("exogenous_prefix is required for method=\"spike-in\"")
            })?;
            bamnado::normalization_factors::spike_in_scale_factors(
                &bam_stats,
                &sample_names,
                &prefix,
            )
            .map_err(anyhow_to_pyerr)?
        } else {
            let filters: Vec<BamReadFilter> = bam_stats
                .iter()
                .map(|stats| build_bam_read_filter(read_filter.clone(), stats))
                .collect::<PyResult<Vec<_>>>()?;

            let counts = bamnado::bin_counts::count_bins(
                &bam_paths,
                &sample_names,
                bin_size,
                &filters,
                use_fragment,
                ignore_scaffold_chromosomes,
            )
            .map_err(anyhow_to_pyerr)?;

            let params = bamnado::normalization_factors::NormFactorParams {
                reference_sample,
                logratio_trim,
                sum_trim,
                exclude_top_percent,
                ..Default::default()
            };
            let scale_factors = bamnado::normalization_factors::compute_scale_factors(
                &counts,
                &method_enum,
                &params,
            )
            .map_err(anyhow_to_pyerr)?;

            if return_counts {
                let bin_chroms = counts.bins.iter().map(|b| b.chrom.clone()).collect();
                let bin_starts = counts.bins.iter().map(|b| b.start).collect();
                let bin_ends = counts.bins.iter().map(|b| b.end).collect();
                counts_data = Some((counts.counts, bin_chroms, bin_starts, bin_ends));
            }

            scale_factors
        };

        let (counts, bin_chroms, bin_starts, bin_ends) = match counts_data {
            Some((c, bc, bs, be)) => (Some(c), bc, bs, be),
            None => (None, Vec::new(), Vec::new(), Vec::new()),
        };

        Ok(ScaleFactorResult {
            method: result.method,
            sample_names: result.sample_names,
            library_sizes: result.library_sizes,
            norm_factors: result.norm_factors,
            scale_factors: result.scale_factors,
            reference_sample: result.reference_sample,
            n_bins_total: result.n_bins_total,
            n_bins_used: result.n_bins_used,
            counts,
            bin_chroms,
            bin_starts,
            bin_ends,
        })
    }

    #[pyfunction]
    fn __version__() -> &'static str {
        env!("CARGO_PKG_VERSION")
    }

    /// Run the bamnado CLI with `argv` and return its exit code.
    #[pyfunction]
    fn _cli_main(py: Python<'_>, argv: Vec<String>) -> i32 {
        py.detach(|| bamnado::cli::run_cli(argv))
    }
}
