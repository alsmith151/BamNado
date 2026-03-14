# BamNado

High-performance tools and utilities for working with BAM and BigWig files in modern genomics workflows. BamNado is written in Rust for speed and low memory use and provides both a command-line interface and Python bindings.

---

## Overview

BamNado is designed for efficient, streaming manipulation of BAM files and signal tracks. It focuses on fast coverage generation, flexible filtering, and lightweight post-processing of bedGraph and BigWig data.

Common use cases include:

- Rapid generation of coverage tracks from large BAM files
- Filtering reads by tags or barcodes to produce targeted BigWigs
- Fragment-aware coverage for ATAC-seq and related assays
- BigWig comparison and aggregation across samples
- Post-processing of binned signal tracks for visualization

BamNado is useful in a range of workflows, including single-cell and Micro-Capture-C (MCC), but is not limited to those applications.

---

## Features

- High-performance, streaming implementations in Rust
- Cross-platform support (Linux, macOS, Windows)
- BAM → bedGraph / BigWig coverage generation
- Fragment-aware and strand-specific pileups
- Read filtering by mapping quality, length, strand, fragment size, tags, and barcodes
- BigWig comparison (subtraction, ratio, log-ratio)
- BigWig aggregation (sum, mean, median, min, max)
- `collapse-bedgraph` utility to merge adjacent bins with identical scores
- Python bindings for selected functionality

---

## Installation

### Pre-built binaries (recommended)

Download the appropriate binary from the [releases page](https://github.com/alsmith151/BamNado/releases).

After downloading:

```bash
chmod +x bamnado
./bamnado --version
```

(Optional) install system-wide:

```bash
sudo cp bamnado /usr/local/bin/
```

---

### Docker

```bash
docker pull ghcr.io/alsmith151/bamnado:latest
docker run --rm ghcr.io/alsmith151/bamnado:latest --help
```

Images are available for `linux/amd64` and `linux/arm64`.

---

### Cargo

If you have Rust installed:

```bash
cargo install bamnado
```

---

### Build from source

```bash
git clone https://github.com/alsmith151/BamNado.git
cd BamNado
cargo build --release
```

---

### Optional dependency: samtools

`samtools` is not required but is **strongly recommended** if your BAM files have non-standard or incomplete headers (e.g. files produced by CellRanger). BamNado automatically falls back to `samtools view -H` to parse the header when the built-in parser fails. Without `samtools` on your `PATH`, BamNado will error on such files.

Install via conda or your system package manager:

```bash
conda install -c bioconda samtools
# or
brew install samtools
```

---

## Python Interface

BamNado provides Python bindings for selected high-performance operations and is available directly from PyPI.

### Python Installation

```bash
pip install bamnado
# or
uv pip install bamnado
```

### ReadFilter

All read filtering options are controlled through the `ReadFilter` class:

| Parameter | Type | Default | Description |
| --------- | ---- | ------- | ----------- |
| `min_mapq` | `int` | `0` | Minimum mapping quality score |
| `proper_pair` | `bool` | `True` | Keep only properly paired reads |
| `min_length` | `int` | `0` | Minimum read length (bp) |
| `max_length` | `int` | `1000` | Maximum read length (bp) |
| `strand` | `str` | `"both"` | Strand to keep: `"forward"` / `"fwd"` / `"+"`, `"reverse"` / `"rev"` / `"-"`, or `"both"` |
| `min_fragment_length` | `int \| None` | `None` | Minimum insert size / TLEN (bp); requires paired-end data |
| `max_fragment_length` | `int \| None` | `None` | Maximum insert size / TLEN (bp); requires paired-end data |
| `blacklist_bed` | `str \| None` | `None` | Path to a BED file of regions to exclude |
| `whitelisted_barcodes` | `list[str] \| None` | `None` | Cell barcodes (CB tag) to include |
| `read_group` | `str \| None` | `None` | Read group (RG tag) to keep |
| `filter_tag` | `str \| None` | `None` | Two-character SAM tag to filter on (e.g. `"VP"`) |
| `filter_tag_value` | `str \| None` | `None` | Required string value for `filter_tag` |

A `ValueError` is raised if `min_fragment_length` or `max_fragment_length` is set on a single-end BAM file.

### Example

```python
import bamnado
import numpy as np

# Basic coverage — default filter settings
signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=False,
    ignore_scaffold_chromosomes=True,
)
print(f"Mean coverage: {np.mean(signal):.3f}")

# Forward-strand nucleosome-free region coverage (100–200 bp fragments)
nfr_filter = bamnado.ReadFilter(
    strand="forward",
    min_fragment_length=100,
    max_fragment_length=200,
    min_mapq=20,
)
nfr_signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=10,
    scale_factor=1.0,
    use_fragment=True,
    ignore_scaffold_chromosomes=True,
    read_filter=nfr_filter,
)

# Tag-filtered coverage (e.g. MCC viewpoint)
vp_filter = bamnado.ReadFilter(
    filter_tag="VP",
    filter_tag_value="BCL2",
    min_mapq=30,
)
vp_signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=True,
    ignore_scaffold_chromosomes=True,
    read_filter=vp_filter,
)
```

---

## Command-line usage

List available commands:

```bash
bamnado --help
```

Get help for a specific command:

```bash
bamnado <command> --help
```

### Available commands

- `bam-coverage` – generate coverage from a BAM file
- `multi-bam-coverage` – coverage from multiple BAMs
- `split` – split BAMs based on filters (e.g. barcodes)
- `split-exogenous` – split endogenous vs exogenous reads
- `modify` – apply transformations and filters to BAMs
- `bigwig-compare` – compare two BigWigs
- `bigwig-aggregate` – aggregate multiple BigWigs
- `collapse-bedgraph` – merge adjacent bedGraph bins with identical scores

---

## Read filtering

All coverage commands share a common set of read filter flags:

| Flag | Default | Description |
| ---- | ------- | ----------- |
| `--strand` | `both` | Include only `forward`, `reverse`, or `both` strands |
| `--proper-pair` | off | Keep only properly-paired reads |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--min-length` | 20 | Minimum read sequence length (bp) |
| `--max-length` | 1000 | Maximum read sequence length (bp) |
| `--min-fragment-length` | — | Minimum insert size / TLEN (bp); requires paired-end data |
| `--max-fragment-length` | — | Maximum insert size / TLEN (bp); requires paired-end data |
| `--blacklisted-locations` | — | BED file of regions to exclude |
| `--whitelisted-barcodes` | — | Text file of cell barcodes to keep (one per line) |
| `--read-group` | — | Keep only reads belonging to this read group |
| `--filter-tag` / `--filter-tag-value` | — | Keep reads where SAM tag equals the given value |

Fragment length filtering operates on the SAM `TLEN` field and is only meaningful for paired-end BAMs. BamNado will return an error if these flags are used with a single-end file.

---

## Example: BAM coverage

```bash
bamnado bam-coverage \
  --bam input.bam \
  --output output.bedgraph \
  --bin-size 100 \
  --norm-method rpkm \
  --scale-factor 1.5 \
  --use-fragment \
  --proper-pair \
  --min-mapq 30
```

---

## Example: strand- and fragment-length-filtered coverage

Useful for isolating nucleosome-free regions in ATAC-seq data:

```bash
bamnado bam-coverage \
  --bam atac.bam \
  --output nfr_forward.bw \
  --bin-size 10 \
  --use-fragment \
  --strand forward \
  --min-fragment-length 100 \
  --max-fragment-length 200 \
  --min-mapq 20
```

---

## Example: tag-filtered BigWig generation

```bash
bamnado bam-coverage \
  --bam input.bam \
  --output BCL2.bw \
  --bin-size 50 \
  --filter-tag "VP" \
  --filter-tag-value "BCL2" \
  --use-fragment \
  --min-mapq 30
```

---

## BigWig comparison

```bash
bamnado bigwig-compare \
  --bw1 sample1.bw \
  --bw2 sample2.bw \
  --comparison log-ratio \
  --pseudocount 1e-3 \
  -o output.bw
```

---

## BigWig aggregation

```bash
bamnado bigwig-aggregate \
  --bigwigs sample1.bw sample2.bw sample3.bw \
  --method mean \
  -o aggregated.bw
```

---

## collapse-bedgraph

```bash
bamnado collapse-bedgraph \
  --input signal.bedgraph \
  --output signal.collapsed.bedgraph
```

---

## Development

```bash
cargo build --release
cargo test
```

---

## License

Apache-2.0 OR MIT
