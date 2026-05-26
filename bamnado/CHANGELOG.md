# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.1](https://github.com/alsmith151/BamNado/compare/v0.8.0...v0.8.1) (2026-05-26)


### Bug Fixes

* resolve 23 correctness and safety bugs across core modules ([#90](https://github.com/alsmith151/BamNado/issues/90)) ([e06a4ee](https://github.com/alsmith151/BamNado/commit/e06a4eebd48bf191646ad163827b7509e0fdb864))

## [0.8.0](https://github.com/alsmith151/BamNado/compare/v0.7.0...v0.8.0) (2026-05-26)


### Features

* add JSON output options to infer-scale subcommand ([8927b58](https://github.com/alsmith151/BamNado/commit/8927b58a8c127f12cb336c435265c112e30ec0c2))

## [0.7.0](https://github.com/alsmith151/BamNado/compare/v0.6.1...v0.7.0) (2026-05-20)


### Features

* add scale factors to BigWig comparison and aggregation options ([#81](https://github.com/alsmith151/BamNado/issues/81)) ([fa53985](https://github.com/alsmith151/BamNado/commit/fa53985f19fa42efe10fc8dd43689f43caf67dce))

## [Unreleased]

### Added

- Comprehensive documentation (docstrings) for all structs and methods in `bamnado/src/`.
- New Python interface for BamNado functionality, allowing usage of Rust-optimized tools directly from Python.
- New `compare-bigwigs` CLI command to compare two BigWig files using `subtraction`, `ratio`, or `log-ratio`.

### Fixed

- Corrected placement of docstrings in `bam_utils.rs` and `read_filter.rs` to ensure `cargo doc` generates correct documentation.

### Changed

- Refined CLI help output with clearer command descriptions, grouped option sections, and command examples.
- Cleaned up several long option names while keeping backwards-compatible aliases, including `--normalize`, `--fragment-counts`, `--blacklist`, `--barcode-allowlist`, `--tag`, and fragment-length flags.

## [0.3.1] - 2025-07-09

### Added

- Initial release with BAM file manipulation tools
- Support for single cell and MCC use cases
- Cross-platform binary builds (Linux, macOS, Windows)
- Rust-based implementation for high performance
- Pure Rust workflow with automatic version bumping
- GitHub Actions CI/CD pipeline

[Unreleased]: https://github.com/alsmith151/BamNado/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/alsmith151/BamNado/releases/tag/v0.3.1
