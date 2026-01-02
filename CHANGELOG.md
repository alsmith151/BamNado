# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Comprehensive documentation (docstrings) for all structs and methods in `bamnado/src/`.
- New Python interface for BamNado functionality, allowing usage of Rust-optimized tools directly from Python.
- New `compare-bigwigs` CLI command to compare two BigWig files using `subtraction`, `ratio`, or `log-ratio`.

### Fixed

- Corrected placement of docstrings in `bam_utils.rs` and `read_filter.rs` to ensure `cargo doc` generates correct documentation.

### Changed

- Changes in existing functionality

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
