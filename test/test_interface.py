import pytest
import os
import numpy as np
import shutil
import subprocess
from pathlib import Path
from bamnado import get_signal_for_chromosome, compute_scale_factors
from bamnado.cli import main

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


def test_compute_scale_factors_tmm(bam_file):
    """Test computing TMM scale factors for two (identical) samples."""
    result = compute_scale_factors([bam_file, bam_file], method="tmm")

    assert len(result.sample_names) == 2
    assert len(result.scale_factors) == 2
    assert len(result.norm_factors) == 2
    assert len(result.library_sizes) == 2
    for factor in result.scale_factors:
        assert factor == pytest.approx(1.0)


def test_compute_scale_factors_return_counts(bam_file):
    """Test that return_counts=True populates a 2D bin count matrix."""
    result = compute_scale_factors(
        [bam_file, bam_file], method="tmm", return_counts=True
    )

    counts = result.counts()
    assert counts is not None
    assert counts.ndim == 2
    assert counts.dtype == np.uint64
    assert counts.shape[1] == 2
    assert counts.shape[0] == len(result.bin_chroms)
    assert counts.shape[0] == len(result.bin_starts)
    assert counts.shape[0] == len(result.bin_ends)


def test_compute_scale_factors_no_counts_by_default(bam_file):
    """Test that counts() is None unless return_counts=True."""
    result = compute_scale_factors([bam_file, bam_file], method="cpm")
    assert result.counts() is None


def test_compute_scale_factors_invalid_method(bam_file):
    """Test that an invalid method string raises ValueError."""
    with pytest.raises(ValueError):
        compute_scale_factors([bam_file, bam_file], method="not-a-method")


def test_compute_scale_factors_nonexistent_bam():
    """Test that a nonexistent BAM file raises RuntimeError."""
    with pytest.raises(RuntimeError):
        compute_scale_factors(["does_not_exist.bam"], method="cpm")


def test_cli_help_and_parse_errors_do_not_exit_interpreter():
    assert main(["--help"]) == 0
    assert main(["nonesuch"]) == 2


def test_cli_can_run_twice_and_restore_sigint_handler():
    import signal

    previous = signal.getsignal(signal.SIGINT)
    assert main(["--version"]) == 0
    assert main(["--version"]) == 0
    assert signal.getsignal(signal.SIGINT) is previous


@pytest.mark.parametrize("args", [["--version"], ["bam-coverage", "--help"], ["nonesuch"]])
def test_installed_python_cli_matches_rust_binary(args):
    """The installed console script must be an exact second CLI front door."""
    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / "target" / "debug" / "bamnado"
    python_cli = shutil.which("bamnado")
    if not binary.exists() or python_cli is None:
        pytest.skip("requires the Rust binary and installed bamnado console script")

    rust_result = subprocess.run(
        [str(binary), *args], capture_output=True, text=True, check=False
    )
    python_result = subprocess.run(
        [python_cli, *args], capture_output=True, text=True, check=False
    )

    assert python_result.returncode == rust_result.returncode
    assert python_result.stdout == rust_result.stdout
    assert python_result.stderr == rust_result.stderr
