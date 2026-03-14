# BamNado

High-performance BAM file processing library for genomics, with a Rust core and Python bindings.

## Workspace Layout

```
BamNado/
├── bamnado/          # Core Rust library + CLI binary
│   └── src/
│       ├── lib.rs                  # Library entry point / public API
│       ├── main.rs                 # CLI (clap-based)
│       ├── coverage_analysis.rs    # BamPileup — parallel pileup/coverage engine
│       ├── read_filter.rs          # BamReadFilter — per-read filtering logic
│       ├── genomic_intervals.rs    # IntervalMaker — BAM record → genomic interval
│       ├── bam_utils.rs            # BamStats, helpers, Iv type alias
│       └── signal_normalization.rs # Raw / CPM / RPKM normalisation
├── bamnado-python/   # PyO3 Python extension (cdylib)
│   ├── src/lib.rs    # Python bindings — exposes get_signal_for_chromosome()
│   └── python/bamnado/__init__.py
├── Cargo.toml        # Workspace root
└── pyproject.toml    # maturin build config for the Python package
```

## Key Types

| Type | Location | Purpose |
|------|----------|---------|
| `BamPileup` | `coverage_analysis.rs` | Parallel per-chromosome coverage computation |
| `BamReadFilter` | `read_filter.rs` | Multi-criterion read filter (strand, MAPQ, length, fragment length, blacklist, barcode…) |
| `IntervalMaker` | `genomic_intervals.rs` | Converts BAM records to `Iv` intervals (fragment or read mode, with optional shift/truncate) |
| `Iv` | `bam_utils.rs` | `Interval<usize, u32>` from `rust-lapper` |
| `Shift` / `Truncate` | `genomic_intervals.rs` | 5′/3′ coordinate adjustments (e.g. Tn5) |

## Python API

```python
from bamnado import ReadFilter, get_signal_for_chromosome

# All filter options are bundled in ReadFilter
rf = ReadFilter(
    min_mapq=20,
    strand="forward",           # "forward"/"fwd"/"+", "reverse"/"rev"/"-", "both" (default)
    min_fragment_length=100,    # insert-size filter (bp); requires paired-end data
    max_fragment_length=200,
    blacklist_bed="blacklist.bed",
    whitelisted_barcodes=["ACGT-1", "TTGA-1"],
    read_group="RG1",
    filter_tag="VP",
    filter_tag_value="BCL2",
)

signal = get_signal_for_chromosome(
    bam_path="sample.bam",
    chromosome_name="chr1",
    bin_size=50,
    scale_factor=1.0,
    use_fragment=True,
    ignore_scaffold_chromosomes=True,
    read_filter=rf,             # optional; None uses BamReadFilter::default()
)
# Returns: numpy float32 array, length = chromosome size
```

## CLI Filter Flags

All subcommands that accept `FilterOptions` support:

| Flag | Default | Description |
|------|---------|-------------|
| `--strand` | `both` | `forward`, `reverse`, or `both` |
| `--proper-pair` | off | Keep only properly-paired reads |
| `--min-mapq` | 20 | Minimum mapping quality |
| `--min-length` | 20 | Minimum read sequence length (bp) |
| `--max-length` | 1000 | Maximum read sequence length (bp) |
| `--min-fragment-length` | — | Minimum insert size / TLEN (bp) |
| `--max-fragment-length` | — | Maximum insert size / TLEN (bp) |
| `--blacklisted-locations` | — | BED file of regions to exclude |
| `--whitelisted-barcodes` | — | Text file of cell barcodes (one per line) |
| `--read-group` | — | Keep only this RG tag value |
| `--filter-tag` / `--filter-tag-value` | — | Keep reads where TAG == VALUE |

Example — nucleosome-free region pileup (insert size 100–200 bp, forward strand only):
```
bamnado coverage --bam sample.bam --strand forward \
    --min-fragment-length 100 --max-fragment-length 200 \
    --use-fragment --bin-size 10
```

## BamReadFilter parameter order

```rust
BamReadFilter::new(
    strand,              // bio_types::strand::Strand
    proper_pair,         // bool
    min_mapq,            // Option<u8>
    min_length,          // Option<u32>  — read sequence length
    max_length,          // Option<u32>  — read sequence length
    read_group,          // Option<String>
    blacklisted_locs,    // Option<HashMap<usize, Lapper<usize, u32>>>
    whitelisted_barcodes,// Option<HashSet<String>>
    filter_tag,          // Option<String>
    filter_tag_value,    // Option<String>
    min_fragment_length, // Option<u32>  — template/insert size (TLEN)
    max_fragment_length, // Option<u32>  — template/insert size (TLEN)
)
```

## Build

```bash
# Rust library + CLI
cargo build --release

# Python wheel (requires maturin)
cd bamnado-python
maturin develop          # editable install
maturin build --release  # wheel
```

## Notes

- BAM files must be indexed (`.bai`).
- Pileup is parallelised with Rayon across genomic chunks; chunk size is auto-tuned from `BamStats`.
- Fragment length filtering uses the SAM `TLEN` field; it is only meaningful for paired-end data. Filtering is applied before the fragment interval is computed, so it works correctly in both read and fragment coverage modes.
- Strand filtering uses the SAM reverse-complement flag; applies to both fragment and read modes.
