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
- Read filtering by mapping quality, length, tags, and barcodes
- BigWig comparison (subtraction, ratio, log-ratio)
- BigWig aggregation (sum, mean, median, min, max)
- `collapse-bedgraph` utility to merge adjacent bins with identical scores
- Python bindings for selected functionality

---

## Installation

### Pre-built binaries (recommended)

Download the appropriate binary from:
https://github.com/alsmith151/BamNado/releases

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

## Python Interface

BamNado provides Python bindings for selected high-performance operations and is available directly from PyPI.

### Installation

```bash
pip install bamnado
# or
uv pip install bamnado
```

### Example

```python
import bamnado
import numpy as np

signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=False,
    ignore_scaffold_chromosomes=True
)

print(f"Mean coverage: {np.mean(signal)}")
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
