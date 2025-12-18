import numpy as np

def get_signal_for_chromosome(
    bam_path: str,
    chromosome_name: str,
    bin_size: int,
    scale_factor: float,
    use_fragment: bool,
    ignore_scaffold_chromosomes: bool,
) -> np.ndarray:
    """Get signal for a specific chromosome from a BAM file.
    Args:
        bam_path (str): Path to the input BAM file.
        chromosome_name (str): Name of the chromosome to process.
        bin_size (int): Size of the bins for signal calculation.
        scale_factor (float): Factor to scale the signal.
        use_fragment (bool): Whether to use fragment length for signal calculation.
        ignore_scaffold_chromosomes (bool): Whether to ignore scaffold chromosomes.
    Returns:
        np.ndarray: Array representing the signal for the specified chromosome.
    """
    ...

def __version__() -> str:
    """Get the version of the bamnado package."""
    ...
