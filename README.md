# Installation and Usage of the `bamnado` Package

## Installation

To install the `bamnado` package, follow these steps:

1. **Install Rust**: Ensure you have Rust installed on your system. You can install it using `rustup`:

    ```bash
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
    ```

2. **Clone the Repository**: Clone the `bamnado` repository to your local machine:

    ```bash
    git clone https://github.com/your-repo/bamnado.git
    cd bamnado
    ```


3. **Install**: Ensure you have Python and `maturin` installed:

    ```bash
    pip install maturin
    ```

    Then, build and install the Python package:

    ```bash
    pip install .
    ```

## Usage

The `bamnado` package provides a command-line interface (CLI) for processing BAM files. Below is an example of its usage.

### Example: Filtering Reads from a BAM File

#### Command

The following command filters reads from a BAM file based on specific criteria and writes the coverage to a `bedGraph` file:

```bash
bamnado bam-coverage \
  --bam input.bam \
  --output output.bedgraph \
  --bin-size 100 \
  --norm-method rpkm \
  --scale-factor 1.5 \
  --use-fragment \
  --proper-pair \
  --min-mapq 30 \
  --min-length 50 \
  --max-length 500 \
  --blacklisted-locations blacklist.bed \
  --whitelisted-barcodes barcodes.txt
```

#### Explanation of Options

- `--bam`: Path to the input BAM file.
- `--output`: Path to the output file (e.g., `bedGraph` or `BigWig`).
- `--bin-size`: Size of genomic bins for coverage calculation.
- `--norm-method`: Normalization method (`raw`, `rpkm`, or `cpm`).
- `--scale-factor`: Scaling factor for normalization.
- `--use-fragment`: Use fragments instead of individual reads for counting.
- `--proper-pair`: Include only properly paired reads.
- `--min-mapq`: Minimum mapping quality for reads to be included.
- `--min-length`: Minimum read length.
- `--max-length`: Maximum read length.
- `--blacklisted-locations`: Path to a BED file specifying regions to exclude.
- `--whitelisted-barcodes`: Path to a file with barcodes to include.

#### Output

The output file (`output.bedgraph`) will contain the normalized coverage data for the BAM file, filtered based on the specified criteria.

### Additional Commands

#### Multi-BAM Coverage

To calculate coverage from multiple BAM files:

```bash
bamnado multi-bam-coverage \
  --bams file1.bam file2.bam \
  --output output.bedgraph \
  --bin-size 100 \
  --norm-method rpkm \
  --scale-factor 1.5 \
  --use-fragment \
  --proper-pair \
  --min-mapq 30 \
  --min-length 50 \
  --max-length 500
```

#### Split BAM File into Endogenous and Exogenous Reads

To split a BAM file into endogenous and exogenous reads:

```bash
bamnado split-exogenous \
  --input input.bam \
  --output output_prefix \
  --exogenous-prefix "exo_" \
  --stats stats.json \
  --allow-unknown-mapq \
  --proper-pair \
  --min-mapq 30 \
  --min-length 50 \
  --max-length 500
```

#### Split BAM File by Cell Barcodes

To split a BAM file based on cell barcodes:

```bash
bamnado split \
  --input input.bam \
  --output output_prefix \
  --whitelisted-barcodes barcodes.txt \
  --proper-pair \
  --min-mapq 30 \
  --min-length 50 \
  --max-length 500
```

## Help

For more details on available commands and options, run:

```bash
bamnado --help
```