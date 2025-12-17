# BamNado

High-performance tools and utilities for manipulation of BAM files for specialized use cases, including single cell and MCC (Multi-modal cellular characterization) workflows.

## Overview

BamNado is a Rust-based toolkit designed to handle complex BAM file operations that are common in modern genomics workflows, particularly in single-cell and multi-modal cellular characterization experiments. It provides efficient, cross-platform tools for coverage calculation, read filtering, file splitting, and various BAM file transformations.

## Python Interface

BamNado also provides a Python interface for direct access to its high-performance BAM processing capabilities.

### Installation

You can install the Python package directly from the source using `pip` or `uv`:

```bash
pip install .
# or
uv pip install .
```

### Usage

```python
import bamnado
import numpy as np

# Get coverage signal for a chromosome
signal = bamnado.get_signal_for_chromosome(
    bam_path="path/to/file.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=False,
    ignore_scaffold_chromosomes=True
)

# signal is a numpy array of floats
print(f"Mean coverage: {np.mean(signal)}")
```

## Installation

BamNado can be installed in several ways. Choose the method that best fits your needs:

### Method 1: Pre-built Binaries (Recommended)

The easiest way to get started is to download a pre-compiled binary from our [releases page](https://github.com/alsmith151/BamNado/releases).

#### Available Platforms

| Platform | Architecture | File Name |
|----------|-------------|-----------|
| Linux | x86_64 | `bamnado-x86_64-unknown-linux-gnu.tar.gz` |
| macOS | Intel (x86_64) | `bamnado-x86_64-apple-darwin.tar.gz` |
| macOS | Apple Silicon (ARM64) | `bamnado-aarch64-apple-darwin.tar.gz` |
| Windows | x86_64 | `bamnado-x86_64-pc-windows-msvc.zip` |

#### Installation Steps

1. **Download the binary**

   Go to the [releases page](https://github.com/alsmith151/BamNado/releases) and download the appropriate file for your system.

2. **Extract the archive**

   **Linux/macOS:**

   ```bash
   tar -xzf bamnado-*.tar.gz
   ```

   **Windows:**
   - Right-click the zip file and select "Extract All"
   - Or use your preferred extraction tool (7-Zip, WinRAR, etc.)

3. **Make executable** (Linux/macOS only)

   ```bash
   chmod +x bamnado
   ```

4. **Test the installation**

   ```bash
   ./bamnado --version
   ```

   You should see output like: `bamnado 0.3.2`

5. **Install system-wide** (optional but recommended)

   **Option A: System-wide installation (requires admin privileges)**

   ```bash
   # Linux/macOS
   sudo cp bamnado /usr/local/bin/

   # Windows (as Administrator)
   # Copy bamnado.exe to C:\Windows\System32\ or add to PATH
   ```

   **Option B: User-local installation (no admin required)**

   ```bash
   # Linux/macOS
   mkdir -p ~/.local/bin
   cp bamnado ~/.local/bin/

   # Add to your shell profile if not already in PATH
   echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
   # or for zsh users:
   echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc

   # Reload your shell or run:
   source ~/.bashrc  # or ~/.zshrc
   ```

6. **Verify system installation**

   Open a new terminal and run:

   ```bash
   bamnado --version
   ```

#### Troubleshooting Pre-built Binaries

##### Linux: "No such file or directory" error

- Your system might be missing required libraries. Try:

  ```bash
  ldd bamnado  # Check dependencies
  ```

- For older Linux distributions, you may need to build from source.

##### macOS: "Cannot be opened because the developer cannot be verified"

- Run: `xattr -d com.apple.quarantine bamnado`
- Or go to System Preferences → Security & Privacy and allow the app

##### Windows: "Windows protected your PC"

- Click "More info" → "Run anyway"
- Or add an exception in Windows Defender

### Method 2: Install via Cargo

If you have Rust and Cargo installed, you can install BamNado directly from crates.io:

```bash
cargo install bamnado
```

**Prerequisites:**

- Rust 1.70+ (install from [rustup.rs](https://rustup.rs/))
- Cargo (comes with Rust)

**Advantages:**

- Always gets the latest published version
- Automatically handles dependencies
- Works on any platform supported by Rust

### Method 3: Build from Source

For the latest development version or if pre-built binaries don't work on your system:

#### Prerequisites

- Rust 2024 edition or later
- Git
- C compiler (for some dependencies)

**Install Rust if you haven't already:**

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

#### Build Steps

1. **Clone the repository**

   ```bash
   git clone https://github.com/alsmith151/BamNado.git
   cd BamNado
   ```

2. **Build the project**

   ```bash
   # Debug build (faster compilation, slower execution)
   cargo build

   # Release build (slower compilation, faster execution - recommended)
   cargo build --release
   ```

3. **Test the build**

   ```bash
   # For debug build
   ./target/debug/bamnado --version

   # For release build
   ./target/release/bamnado --version
   ```

4. **Install system-wide** (optional)

   ```bash
   # Install from source
   cargo install --path .

   # Or manually copy the binary
   sudo cp target/release/bamnado /usr/local/bin/
   ```

#### Build Troubleshooting

##### Common Issues

##### Error: "linker 'cc' not found"

- **Ubuntu/Debian:** `sudo apt install build-essential`
- **CentOS/RHEL:** `sudo yum groupinstall "Development Tools"`
- **macOS:** Install Xcode Command Line Tools: `xcode-select --install`
- **Windows:** Install Visual Studio Build Tools or use WSL

##### Error: "failed to run custom build command for 'openssl-sys'"

- **Ubuntu/Debian:** `sudo apt install libssl-dev pkg-config`
- **CentOS/RHEL:** `sudo yum install openssl-devel pkgconf-pkg-config`
- **macOS:** Usually works out of the box with Homebrew
- **Windows:** Consider using the pre-built binaries instead

### Quick Start Verification

After installation, verify everything works:

```bash
# Check version
bamnado --version

# See available commands
bamnado --help

# Test with a simple command (replace with your BAM file)
bamnado bam-coverage --bam /path/to/your/file.bam --output test.bedgraph
```

## Usage

### Available Commands

BamNado provides several commands for different BAM file operations:

- `bam-coverage` - Calculate coverage from a BAM file and write to a bedGraph or bigWig file
- `multi-bam-coverage` - Calculate coverage from multiple BAM files and write to a bedGraph or bigWig file
- `split-exogenous` - Split a BAM file into endogenous and exogenous reads
- `split` - Split a BAM file based on a set of defined filters
- `modify` - Modify BAM files with various transformations

For detailed help on any command, use:

```bash
bamnado <command> --help
```

### Example: Calculating Coverage from a BAM File

#### Command

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
- `--min-mapq`: Minimum mapping quality for reads to be included (default: 20).
- `--min-length`: Minimum read length (default: 20).
- `--max-length`: Maximum read length (default: 1000).
- `--blacklisted-locations`: Path to a BED file specifying regions to exclude.
- `--whitelisted-barcodes`: Path to a file with barcodes to include.
- `--strand`: Filter reads based on strand (both, forward, reverse).
- `--shift`: Shift options for the pileup (default: 0,0,0,0).
- `--truncate`: Truncate options for the pileup.
- `--ignore-scaffold`: Ignore scaffold chromosomes.
- `--read-group`: Selected read group.

#### Output

The output file (`output.bedgraph`) will contain the normalized coverage data for the BAM file, filtered based on the specified criteria. BigWig files can also be generated by specifying the `--output` option with a `.bw` extension.

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

#### Modify BAM Files

To modify BAM files with various transformations:

```bash
bamnado modify \
  --input input.bam \
  --output output_prefix \
  --proper-pair \
  --min-mapq 30 \
  --min-length 50 \
  --max-length 500 \
  --tn5-shift
```

The `modify` command supports various filtering options and transformations like Tn5 shifting for ATAC-seq data processing.

## Help

For more details on available commands and options, run:

```bash
bamnado --help
```

Or for specific command help:

```bash
bamnado <command> --help
```

## Features

- **High Performance**: Built in Rust for maximum speed and memory efficiency
- **Cross-platform**: Available for Linux, macOS, and Windows
- **Multiple Output Formats**: Support for bedGraph and BigWig output formats
- **Flexible Filtering**: Comprehensive read filtering options including mapping quality, read length, proper pairs, and more
- **Single Cell Support**: Built-in support for cell barcode-based operations
- **MCC Workflows**: Specialized tools for Multi-modal Cellular Characterization
- **Strand-specific Analysis**: Support for strand-specific coverage calculations
- **Blacklist/Whitelist Support**: Region and barcode filtering capabilities

## Development

### Requirements

- Rust 2024 edition or later
- Cargo package manager

### Building from Source

```bash
git clone https://github.com/alsmith151/BamNado.git
cd BamNado
cargo build --release
```

### Running Tests

```bash
cargo test
```

### Pre-commit Hooks

This project uses pre-commit hooks to ensure code quality and consistency. The hooks run the same checks as the CI workflow:

- Code formatting (`cargo fmt`)
- Linting (`cargo clippy`)
- Basic checks (`cargo check`)
- Tests (`cargo test` on push)

#### Quick Setup

Run the setup script to install and configure pre-commit hooks:

```bash
./setup-precommit.sh
```

#### Manual Setup

If you prefer to set up pre-commit manually:

```bash
# Install pre-commit (choose one method)
pip install pre-commit
# or: brew install pre-commit
# or: conda install -c conda-forge pre-commit

# Install the hooks
pre-commit install
pre-commit install --hook-type pre-push

# Test the setup
pre-commit run --all-files
```

#### Configuration Options

Two pre-commit configurations are available:

- `.pre-commit-config.yaml` - Full checks including `cargo check` on every commit
- `.pre-commit-config-fast.yaml` - Faster setup with formatting/linting only, tests on push

To use the fast configuration:

```bash
mv .pre-commit-config.yaml .pre-commit-config-full.yaml
mv .pre-commit-config-fast.yaml .pre-commit-config.yaml
pre-commit install
```

#### Useful Commands

```bash
pre-commit run --all-files       # Run all hooks on all files
pre-commit run cargo-fmt         # Run specific hook
pre-commit autoupdate            # Update hook versions
pre-commit uninstall             # Remove hooks
```

## Release Information

### Version 0.3.1 (2025-07-09)

- Initial public release with comprehensive BAM file manipulation tools
- Support for single cell and MCC (Multi-modal Cellular Characterization) use cases
- Cross-platform binary builds available for Linux, macOS, and Windows
- High-performance Rust implementation
- Complete CI/CD pipeline with automated testing and releases

For detailed changelog information, see [CHANGELOG.md](CHANGELOG.md).

## License

This project is licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
