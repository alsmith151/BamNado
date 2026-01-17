//! # BamNado Library
//!
//! `bamnado` is a high-performance Rust library for manipulating BAM files, designed specifically
//! for specialized genomics workflows such as single-cell ATAC-seq (scATAC-seq) and Multi-modal
//! Cellular Characterization (MCC).
//!
//! ## Core Features
//!
//! *   **Coverage Analysis**: Efficient calculation of coverage signals across chromosomes.
//! *   **BAM Splitting**: Tools for splitting BAM files based on various criteria.
//! *   **Read Filtering**: Advanced filtering capabilities for BAM reads.
//! *   **Signal Normalization**: Methods for normalizing genomic signals.
//! *   **Spike-in Analysis**: Utilities for handling spike-in controls.
//!
//! ## Modules
//!
//! *   [`bam_modifier`]: Tools for modifying BAM records.
//! *   [`bam_splitter`]: Functionality for splitting BAM files.
//! *   [`bam_utils`]: General utilities and statistics for BAM files.
//! *   [`coverage_analysis`]: Core logic for generating coverage tracks.
//! *   [`genomic_intervals`]: Handling of genomic coordinates and intervals.
//! *   [`read_filter`]: Structures and traits for filtering reads.
//! *   [`signal_normalization`]: Normalization strategies for coverage data.
//! *   [`spike_in_analysis`]: Analysis tools for spike-in normalization.

// Rust crates

// Internal modules
pub mod bam_modifier;
pub mod bam_splitter;
pub mod bam_utils;
pub mod bedgraph_utils;
pub mod bigwig_compare;
pub mod coverage_analysis;
pub mod genomic_intervals;
pub mod read_filter;
pub mod signal_normalization;
pub mod spike_in_analysis;

pub use bam_utils::{BamStats, CellBarcodes, CellBarcodesMulti, Strand};
