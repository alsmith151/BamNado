import numpy as np

class ReadFilter:
    """Configuration for filtering BAM reads.

    Controls which reads are included in coverage calculations based on mapping
    quality, read length, strand, paired-end flags, barcodes, and genomic blacklists.
    """

    min_mapq: int
    proper_pair: bool
    min_length: int
    max_length: int
    strand: str
    min_fragment_length: int | None
    max_fragment_length: int | None
    blacklist_bed: str | None
    whitelisted_barcodes: list[str] | None
    read_group: str | None
    filter_tag: str | None
    filter_tag_value: str | None

    def __init__(
        self,
        min_mapq: int = 0,
        proper_pair: bool = False,
        min_length: int = 0,
        max_length: int = 1000,
        strand: str = "both",
        min_fragment_length: int | None = None,
        max_fragment_length: int | None = None,
        blacklist_bed: str | None = None,
        whitelisted_barcodes: list[str] | None = None,
        read_group: str | None = None,
        filter_tag: str | None = None,
        filter_tag_value: str | None = None,
    ) -> None: ...

def get_signal_for_chromosome(
    bam_path: str,
    chromosome_name: str,
    bin_size: int,
    scale_factor: float,
    use_fragment: bool,
    ignore_scaffold_chromosomes: bool,
    read_filter: ReadFilter | None = None,
) -> np.ndarray:
    """Get signal for a specific chromosome from a BAM file.
    Args:
        bam_path (str): Path to the input BAM file.
        chromosome_name (str): Name of the chromosome to process.
        bin_size (int): Size of the bins for signal calculation.
        scale_factor (float): Factor to scale the signal.
        use_fragment (bool): Whether to use fragment length for signal calculation.
        ignore_scaffold_chromosomes (bool): Whether to ignore scaffold chromosomes.
        read_filter (ReadFilter | None): Optional read filter configuration.
    Returns:
        np.ndarray: Array representing the signal for the specified chromosome.
    """
    ...

class ScaleFactorResult:
    """Result of a between-sample scaling factor computation."""

    method: str
    sample_names: list[str]
    library_sizes: list[int]
    norm_factors: list[float]
    scale_factors: list[float]
    reference_sample: str | None
    n_bins_total: int
    n_bins_used: int

    def counts(self) -> np.ndarray | None:
        """The (n_bins, n_samples) bin count matrix, or None if return_counts=False."""
        ...

    @property
    def bin_chroms(self) -> list[str]: ...
    @property
    def bin_starts(self) -> list[int]: ...
    @property
    def bin_ends(self) -> list[int]: ...

def compute_scale_factors(
    bam_paths: list[str],
    method: str = "csaw-background",
    bin_size: int = 10000,
    sample_names: list[str] | None = None,
    exclude_top_percent: float = 5.0,
    reference_sample: str | None = None,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
    exogenous_prefix: str | None = None,
    use_fragment: bool = False,
    ignore_scaffold_chromosomes: bool = True,
    read_filter: ReadFilter | None = None,
    return_counts: bool = False,
) -> ScaleFactorResult:
    """Compute between-sample scaling factors from BAM files.

    Args:
        bam_paths (list[str]): Input BAM file paths, one per sample.
        method (str): Normalization method: "tmm", "csaw-background" (default), "cpm",
            "median-of-ratios", or "spike-in".
        bin_size (int): Genomic bin width in base pairs. Ignored for "spike-in".
        sample_names (list[str] | None): Sample names, same length/order as bam_paths.
            Defaults to each BAM file's stem, disambiguated if duplicated.
        exclude_top_percent (float): Percentage of highest-mean-count bins dropped
            before estimating ("csaw-background" only).
        reference_sample (str | None): Sample name to use as the TMM reference.
        logratio_trim (float): TMM trim fraction for M-values.
        sum_trim (float): TMM trim fraction for A-values.
        exogenous_prefix (str | None): Reference-name prefix identifying spike-in
            sequences. Required for method="spike-in".
        use_fragment (bool): Count fragments (pairs) instead of individual reads.
        ignore_scaffold_chromosomes (bool): Skip scaffold / unplaced chromosomes.
        read_filter (ReadFilter | None): Read filter configuration, applied per sample.
        return_counts (bool): If True, populate ScaleFactorResult.counts() with the
            underlying bin count matrix.

    Returns:
        ScaleFactorResult: The computed scaling factors and (optionally) count matrix.
    """
    ...

def __version__() -> str:
    """Get the version of the bamnado package."""
    ...
