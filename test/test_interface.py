import pytest
import os
import numpy as np
from bamnado import get_signal_for_chromosome

@pytest.fixture
def bam_file():
    test_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(test_dir, "data", "test.bam")

def test_get_signal_for_chromosome(bam_file):
    """Test getting signal for a chromosome."""
    assert os.path.exists(bam_file), f"Test file not found: {bam_file}"

    # Parameters
    # Note: Adjust chromosome name based on what is actually in test.bam
    chromosome = "chr9"
    bin_size = 50
    scale_factor = 1.0
    use_fragment = False
    ignore_scaffold = True

    try:
        signal = get_signal_for_chromosome(
            bam_file,
            chromosome,
            bin_size,
            scale_factor,
            use_fragment,
            ignore_scaffold
        )

        assert isinstance(signal, np.ndarray)
        assert signal.ndim == 1
        assert signal.dtype == np.float32

    except KeyError as e:
        pytest.fail(f"Chromosome {chromosome} not found in {bam_file}: {e}")

def test_get_signal_invalid_chrom(bam_file):
    """Test that requesting a non-existent chromosome raises KeyError."""
    with pytest.raises(KeyError):
        get_signal_for_chromosome(
            bam_file,
            "chr_nonexistent_12345",
            100,
            1.0,
            False,
            True
        )
