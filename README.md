# BamNado

High-performance tools and utilities for working with BAM and BigWig files in modern genomics workflows. BamNado is written in Rust for speed and low memory use and provides both a command-line interface and Python bindings.

## Overview

BamNado is designed for efficient, streaming manipulation of BAM files and signal tracks. It focuses on fast coverage generation, flexible filtering, and lightweight post-processing of bedGraph and BigWig data.

**Common use cases:**

- Rapid generation of coverage tracks from large BAM files
- Filtering reads by tags or barcodes to produce targeted BigWigs
- Fragment-aware coverage for ATAC-seq and related assays
- BigWig comparison and aggregation across samples
- Post-processing of binned signal tracks for visualization

Useful in workflows including single-cell and Micro-Capture-C (MCC), and many others.

**Documentation:** [API docs on docs.rs](https://docs.rs/bamnado/latest/bamnado/)

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

### Docker

```bash
docker pull ghcr.io/alsmith151/bamnado:latest
docker run --rm ghcr.io/alsmith151/bamnado:latest --help
```

Images are available for `linux/amd64` and `linux/arm64`.

### Cargo

If you have Rust installed:

```bash
cargo install bamnado
```


### Build from source

```bash
git clone https://github.com/alsmith151/BamNado.git
cd BamNado
cargo build --release
```

### Optional dependency: samtools

`samtools` is not required but is **strongly recommended** if your BAM files have non-standard or incomplete headers (e.g. files produced by CellRanger). BamNado automatically falls back to `samtools view -H` to parse the header when the built-in parser fails. Without `samtools` on your `PATH`, BamNado will error on such files.

Install via conda or your system package manager:

```bash
conda install -c bioconda samtools
# or
brew install samtools
```

## Python API

BamNado provides Python bindings for high-performance BAM signal generation with flexible read filtering. Install from PyPI:

```bash
pip install bamnado
# or
uv pip install bamnado
```

### Quick start

```python
import bamnado
import numpy as np

# Generate coverage from a BAM file
signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=50,
)
print(f"Mean coverage: {np.mean(signal):.2f}")
```

### ReadFilter class

Customize read filtering with the `ReadFilter` class. All parameters are optional:

```python
import bamnado

# Create a filter for high-quality, properly-paired reads
strict_filter = bamnado.ReadFilter(
    min_mapq=30,
    proper_pair=True,
    min_length=50,
)

# Apply the filter
signal = bamnado.get_signal_for_chromosome(
    bam_path="input.bam",
    chromosome_name="chr1",
    bin_size=100,
    read_filter=strict_filter,
)
```

#### Filter parameters

| Parameter | Type | Default | Description |
| --------- | ---- | ------- | ----------- |
| `min_mapq` | `int` | `0` | Minimum mapping quality |
| `proper_pair` | `bool` | `False` | Require properly paired reads |
| `min_length` | `int` | `0` | Minimum read length (bp) |
| `max_length` | `int` | `1000` | Maximum read length (bp) |
| `strand` | `str` | `"both"` | `"forward"`, `"reverse"`, or `"both"` |
| `min_fragment_length` | `int \| None` | `None` | Minimum insert size (bp); paired-end only |
| `max_fragment_length` | `int \| None` | `None` | Maximum insert size (bp); paired-end only |
| `blacklist_bed` | `str \| None` | `None` | BED file of excluded regions |
| `whitelisted_barcodes` | `list[str] \| None` | `None` | Cell barcodes to include (CB tag) |
| `read_group` | `str \| None` | `None` | Read group to keep (RG tag) |
| `filter_tag` | `str \| None` | `None` | SAM tag to filter on (e.g. `"VP"`) |
| `filter_tag_value` | `str \| None` | `None` | Required value for `filter_tag` |

### Usage examples

#### Basic coverage

```python
import bamnado
import numpy as np

signal = bamnado.get_signal_for_chromosome(
    bam_path="sample.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=False,
)
print(f"Signal shape: {signal.shape}, type: {signal.dtype}")
print(f"Mean: {np.mean(signal):.2f}, Max: {np.max(signal):.2f}")
```

#### Nucleosome-free regions (ATAC-seq)

```python
import bamnado

# Forward-strand reads, 100–200 bp fragments
nfr_filter = bamnado.ReadFilter(
    strand="forward",
    min_fragment_length=100,
    max_fragment_length=200,
    min_mapq=20,
)

nfr_signal = bamnado.get_signal_for_chromosome(
    bam_path="atac.bam",
    chromosome_name="chr1",
    bin_size=10,
    use_fragment=True,
    read_filter=nfr_filter,
)
```

#### Barcode-filtered coverage (single-cell)

```python
import bamnado

# Get coverage for specific cell barcodes
cell_barcodes = ["ACGTACGT-1", "TGCATGCA-1", "AAATAAAA-1"]

barcode_filter = bamnado.ReadFilter(
    whitelisted_barcodes=cell_barcodes,
    min_mapq=25,
)

signal = bamnado.get_signal_for_chromosome(
    bam_path="possorted_genome_bam.bam",
    chromosome_name="chr1",
    bin_size=100,
    use_fragment=True,
    read_filter=barcode_filter,
)
```

#### Tag-filtered coverage (e.g. MCC viewpoint)

```python
import bamnado

# Get coverage for specific viewpoint
vp_filter = bamnado.ReadFilter(
    filter_tag="VP",
    filter_tag_value="BCL2",
    min_mapq=30,
)

vp_signal = bamnado.get_signal_for_chromosome(
    bam_path="mcc.bam",
    chromosome_name="chr1",
    bin_size=50,
    use_fragment=True,
    read_filter=vp_filter,
)
```

#### High-confidence coverage

```python
import bamnado

# High-quality, properly-paired reads
confident_filter = bamnado.ReadFilter(
    min_mapq=30,
    proper_pair=True,
    min_length=50,
    min_fragment_length=100,
    max_fragment_length=500,
)

signal = bamnado.get_signal_for_chromosome(
    bam_path="sample.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1e6,  # CPM normalization
    use_fragment=True,
    read_filter=confident_filter,
)
```

### get_signal_for_chromosome parameters

| Parameter | Type | Default | Description |
| --------- | ---- | ------- | ----------- |
| `bam_path` | `str` | — | Path to indexed BAM file |
| `chromosome_name` | `str` | — | Chromosome to process (e.g. `"chr1"`) |
| `bin_size` | `int` | `50` | Bin width in base pairs |
| `scale_factor` | `float` | `1.0` | Linear scaling factor |
| `use_fragment` | `bool` | `False` | Use fragment (pair) coordinates instead of read |
| `ignore_scaffold_chromosomes` | `bool` | `True` | Skip non-standard chromosomes |
| `read_filter` | `ReadFilter \| None` | `None` | Optional read filter |

**Returns:** `numpy.ndarray` of dtype `float32` with length = chromosome size / bin_size

**Note:** Fragment length filtering (`min_fragment_length`, `max_fragment_length`) requires paired-end BAM files and will raise `ValueError` on single-end data.

## Command-line interface

### Getting help

```bash
bamnado --help                      # List all commands
bamnado <command> --help            # Help for a specific command
```

### Available commands

- **Coverage generation**: `bam-coverage`, `multi-bam-coverage`
- **BAM manipulation**: `split`, `split-exogenous`, `modify`
- **BigWig tools**: `bigwig-compare`, `bigwig-aggregate`, `collapse-bedgraph`

### Read filtering

All coverage commands share common read filter flags:

| Flag | Default | Description |
| ---- | ------- | ----------- |
| `--strand` | `both` | `forward`, `reverse`, or `both` |
| `--proper-pair` | off | Keep only properly-paired reads |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--min-length` | 20 | Minimum read length (bp) |
| `--max-length` | 1000 | Maximum read length (bp) |
| `--min-fragment-length` | — | Minimum insert size (bp); paired-end only |
| `--max-fragment-length` | — | Maximum insert size (bp); paired-end only |
| `--blacklisted-locations` | — | BED file of regions to exclude |
| `--whitelisted-barcodes` | — | Text file of cell barcodes (one per line) |
| `--read-group` | — | Keep only this read group |
| `--filter-tag` / `--filter-tag-value` | — | Keep reads where TAG == VALUE |

### Examples

#### Generate basic coverage

```bash
bamnado bam-coverage \
  --bam sample.bam \
  --output coverage.bedgraph \
  --bin-size 100
```

#### Generate high-quality, normalized coverage

```bash
bamnado bam-coverage \
  --bam sample.bam \
  --output coverage_hq.bw \
  --bin-size 100 \
  --norm-method rpkm \
  --min-mapq 30 \
  --proper-pair
```

#### Extract nucleosome-free regions from ATAC-seq

```bash
bamnado bam-coverage \
  --bam atac.bam \
  --output nfr_forward.bw \
  --bin-size 10 \
  --strand forward \
  --min-fragment-length 100 \
  --max-fragment-length 200 \
  --use-fragment \
  --min-mapq 20
```

#### Filter by cell barcode (single-cell)

```bash
bamnado bam-coverage \
  --bam possorted_genome_bam.bam \
  --output cell_coverage.bw \
  --bin-size 100 \
  --whitelisted-barcodes barcodes.txt \
  --use-fragment
```

#### Filter by SAM tag (MCC viewpoint)

```bash
bamnado bam-coverage \
  --bam mcc.bam \
  --output BCL2_viewpoint.bw \
  --bin-size 50 \
  --filter-tag "VP" \
  --filter-tag-value "BCL2" \
  --use-fragment \
  --min-mapq 30
```

#### Compare two BigWig files

```bash
bamnado bigwig-compare \
  --bw1 sample_treated.bw \
  --bw2 sample_control.bw \
  --comparison log-ratio \
  --pseudocount 1e-3 \
  -o treated_vs_control.bw
```

#### Aggregate multiple BigWig files

```bash
bamnado bigwig-aggregate \
  --bigwigs sample1.bw sample2.bw sample3.bw \
  --method mean \
  -o mean_coverage.bw
```

#### Simplify bedGraph file

```bash
bamnado collapse-bedgraph \
  --input signal.bedgraph \
  --output signal.collapsed.bedgraph
```


## Development

```bash
cargo build --release
cargo test
```


## License

Apache-2.0 OR MIT
