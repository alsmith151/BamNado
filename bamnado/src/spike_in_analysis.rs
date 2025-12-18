//! # Spike-in Analysis Module
//!
//! This module handles the analysis of spike-in controls for normalization purposes.
//! It includes tools for:
//! *   Tracking statistics of endogenous vs. exogenous (spike-in) reads.
//! *   Calculating normalization factors based on spike-in counts.
//! *   Serializing and deserializing spike-in statistics.

use ahash::HashMap;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::path::{Path, PathBuf};
use std::prelude::rust_2021::*;

use anyhow::Result;
use crossbeam::channel::unbounded;
use log::info;
use noodles::{bam, sam};

/// Statistics for the read splitting process.
///
/// Tracks the number of reads assigned to each category (endogenous, exogenous, etc.).
#[derive(Debug, Deserialize)]
pub struct SplitStats {
    filename: String,
    n_total_reads: u64,
    n_unmapped_reads: u64,
    n_qcfail_reads: u64,
    n_duplicate_reads: u64,
    n_secondary_reads: u64,
    n_low_maq: u64,
    n_both_genomes: u64,
    n_exogenous: u64,
    n_endogenous: u64,
}

impl SplitStats {
    /// Create a new `SplitStats` struct.
    fn new(filename: String) -> Self {
        Self {
            filename,
            n_total_reads: 0,
            n_unmapped_reads: 0,
            n_qcfail_reads: 0,
            n_duplicate_reads: 0,
            n_secondary_reads: 0,
            n_low_maq: 0,
            n_both_genomes: 0,
            n_exogenous: 0,
            n_endogenous: 0,
        }
    }

    /// Calculate the normalization factor based on the number of exogenous reads.
    fn spikein_norm_factor(&self) -> f64 {
        // Calculate the ChIP spike-in normalization factor
        // df_counts["scale_factor"] = 1 / (df_counts["spikein_reads"] / 1e6)
        // # if df_counts["norm_factor"] == inf change to 1
        // df_counts["scale_factor"] = df_counts["scale_factor"].replace([np.inf, -np.inf], 1)

        let spikein_reads = self.n_exogenous as f64;
        let scale_factor = 1.0 / (spikein_reads / 1e6);
        if scale_factor.is_infinite() {
            1.0_f64
        } else {
            scale_factor
        }
    }

    /// Calculate the percentage of reads that are exogenous.
    fn spikein_percentage(&self) -> f64 {
        (self.n_exogenous as f64 / self.n_endogenous as f64) * 100.0
    }
}

impl Display for SplitStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Statistics for {}", self.filename)?;
        writeln!(f, "Total reads: {}", self.n_total_reads)?;
        writeln!(f, "Endogenous reads: {}", self.n_endogenous)?;
        writeln!(f, "Exogenous reads: {}", self.n_exogenous)?;
        writeln!(f, "Both genomes reads: {}", self.n_both_genomes)?;
        writeln!(f, "\nFiltered reads:")?;
        writeln!(f, "Unmapped reads: {}", self.n_unmapped_reads)?;
        writeln!(f, "QC fail reads: {}", self.n_qcfail_reads)?;
        writeln!(f, "Duplicate reads: {}", self.n_duplicate_reads)?;
        writeln!(f, "Secondary reads: {}", self.n_secondary_reads)?;
        writeln!(f, "Low mapping quality reads: {}", self.n_low_maq)?;

        Ok(())
    }
}

impl Serialize for SplitStats {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::ser::Serializer,
    {
        let mut state = serializer.serialize_struct("SplitStats", 11)?;
        state.serialize_field("filename", &self.filename)?;
        state.serialize_field("n_total_reads", &self.n_total_reads)?;
        state.serialize_field("n_unmapped_reads", &self.n_unmapped_reads)?;
        state.serialize_field("n_qcfail_reads", &self.n_qcfail_reads)?;
        state.serialize_field("n_duplicate_reads", &self.n_duplicate_reads)?;
        state.serialize_field("n_secondary_reads", &self.n_secondary_reads)?;
        state.serialize_field("n_low_maq", &self.n_low_maq)?;
        state.serialize_field("n_both_genomes", &self.n_both_genomes)?;
        state.serialize_field("n_exogenous", &self.n_exogenous)?;
        state.serialize_field("n_endogenous", &self.n_endogenous)?;
        state.serialize_field("spike_in_factor", &self.spikein_norm_factor())?;
        state.serialize_field("spike_in_percentage", &self.spikein_percentage())?;
        state.end()
    }
}

/// A tool for splitting a BAM file into endogenous and exogenous reads.
///
/// Reads are classified based on the reference sequence they are mapped to.
/// Reference sequences starting with `exogenous_prefix` are considered exogenous.
pub struct BamSplitter {
    // The input BAM file
    input_bam: bam::io::Reader<noodles::bgzf::io::Reader<std::fs::File>>,
    #[allow(dead_code)]
    input_header: sam::Header,

    // The output BAM files
    endogenous_bam: PathBuf,
    exogenous_bam: PathBuf,
    both_bam: PathBuf,
    unmapped_bam: PathBuf,

    // Refseq ids for endogenous and exogenous sequences
    endogenous_refseq_ids: HashMap<usize, String>,
    exogenous_refseq_ids: HashMap<usize, String>,

    // The headers of the output BAM files
    endogenous_header: sam::Header,
    exogenous_header: sam::Header,
    both_header: sam::Header,
    unmapped_header: sam::Header,

    // The statistics for the split
    stats: SplitStats,
}

impl BamSplitter {
    /// Create a new `BamSplitter`.
    ///
    /// # Arguments
    ///
    /// * `input_path` - Path to the input BAM file.
    /// * `output_prefix` - Prefix for the output BAM files.
    /// * `exogenous_prefix` - Prefix for exogenous reference sequences.
    pub fn new(
        input_path: PathBuf,
        output_prefix: PathBuf,
        exogenous_prefix: String,
    ) -> Result<Self> {
        let mut input_bam: bam::io::Reader<noodles::bgzf::io::Reader<std::fs::File>> =
            bam::io::reader::Builder.build_from_path(input_path.clone())?;

        let endogenous_path = output_prefix.with_extension("endogenous.bam");
        let exogenous_path = output_prefix.with_extension("exogenous.bam");
        let both_path = output_prefix.with_extension("both.bam");
        let unmapped_path = output_prefix.with_extension("unmapped.bam");

        let header_input = input_bam.read_header()?;
        let reference_seqs = header_input.reference_sequences().clone();

        // Split reference sequences into endogenous and exogenous based on prefixes present.
        // Endogenous sequences have no prefix, exogenous sequences have a prefix.
        let mut reference_seqs_endogenous = sam::header::ReferenceSequences::new();
        let mut reference_seqs_exogenous = sam::header::ReferenceSequences::new();

        let mut endogenous_refseq_ids = HashMap::default();
        let mut exogenous_refseq_ids = HashMap::default();

        for (i, (name, len)) in reference_seqs.iter().enumerate() {
            let name = name.to_string();
            if name.starts_with(&exogenous_prefix) {
                reference_seqs_exogenous.insert(name.clone().into(), len.clone());
                exogenous_refseq_ids.entry(i).or_insert(name);
            } else {
                reference_seqs_endogenous.insert(name.clone().into(), len.clone());
                endogenous_refseq_ids.entry(i).or_insert(name);
            }
        }

        let header_endogenous = header_input.clone();
        let header_exogenous = header_input.clone();
        let header_both_genomes = header_input.clone();
        let header_unmapped = header_input.clone();

        let stats = SplitStats::new(
            input_path
                .file_name()
                .expect("No file name")
                .to_str()
                .expect("Failed to convert to string")
                .to_string(),
        );

        Ok(Self {
            input_bam,
            input_header: header_input,
            endogenous_bam: endogenous_path,
            exogenous_bam: exogenous_path,
            both_bam: both_path,
            unmapped_bam: unmapped_path,
            endogenous_header: header_endogenous,
            exogenous_header: header_exogenous,
            both_header: header_both_genomes,
            unmapped_header: header_unmapped,
            endogenous_refseq_ids,
            exogenous_refseq_ids,
            stats,
        })
    }

    /// Split the input BAM file into endogenous and exogenous BAM files.
    ///
    /// Reads are processed in parallel and written to the corresponding output files.
    pub fn split(&mut self) -> Result<()> {
        let mut counter = 0;

        let (tx_endogenous, rx_endogenous) = unbounded();
        let (tx_exogenous, rx_exogenous) = unbounded();
        let (tx_both, rx_both) = unbounded();
        let (tx_unmapped, rx_unmapped) = unbounded();

        let endogenous_bam = self.endogenous_bam.clone();
        let exogenous_bam = self.exogenous_bam.clone();
        let both_bam = self.both_bam.clone();
        let unmapped_bam = self.unmapped_bam.clone();

        let endogenous_header = self.endogenous_header.clone();
        let exogenous_header = self.exogenous_header.clone();
        let both_header = self.both_header.clone();
        let unmapped_header = self.unmapped_header.clone();

        // Create writer threads for the output BAM files
        let endogenous_writer = std::thread::spawn(move || {
            let mut writer = bam::io::writer::Builder
                .build_from_path(endogenous_bam)
                .expect("Failed to create writer");
            writer
                .write_header(&endogenous_header)
                .expect("Failed to write header");
            while let Ok(record) = rx_endogenous.recv() {
                writer
                    .write_record(&endogenous_header, &record)
                    .expect("Failed to write record");
            }
        });

        let exogenous_writer = std::thread::spawn(move || {
            let mut writer = bam::io::writer::Builder
                .build_from_path(exogenous_bam)
                .expect("Failed to create writer");
            writer
                .write_header(&exogenous_header)
                .expect("Failed to write header");
            while let Ok(record) = rx_exogenous.recv() {
                writer
                    .write_record(&exogenous_header, &record)
                    .expect("Failed to write record");
            }
        });

        let both_writer = std::thread::spawn(move || {
            let mut writer = bam::io::writer::Builder
                .build_from_path(both_bam)
                .expect("Failed to create writer");
            writer
                .write_header(&both_header)
                .expect("Failed to write header");
            while let Ok(record) = rx_both.recv() {
                writer
                    .write_record(&both_header, &record)
                    .expect("Failed to write record");
            }
        });

        let unmapped_writer = std::thread::spawn(move || {
            let mut writer = bam::io::writer::Builder
                .build_from_path(unmapped_bam)
                .expect("Failed to create writer");
            writer
                .write_header(&unmapped_header)
                .expect("Failed to write header");
            while let Ok(record) = rx_unmapped.recv() {
                writer
                    .write_record(&unmapped_header, &record)
                    .expect("Failed to write record");
            }
        });

        // Read records from the input BAM file and write to the appropriate output BAM file
        let mut record = bam::Record::default();
        while self.input_bam.read_record(&mut record)? > 0 {
            self.stats.n_total_reads += 1;

            counter += 1;
            if counter % (1e6 as usize) == 0 {
                info!("Processed {counter} reads");
            }

            let is_unmapped = record.flags().is_unmapped();
            let is_qcfail = record.flags().is_qc_fail();
            let is_duplicate = record.flags().is_duplicate();
            let is_secondary = record.flags().is_secondary();

            if is_unmapped {
                self.stats.n_unmapped_reads += 1;
                tx_unmapped
                    .send(record.clone())
                    .expect("Failed to send record");
            } else if is_qcfail {
                self.stats.n_qcfail_reads += 1;
                tx_unmapped
                    .send(record.clone())
                    .expect("Failed to send record");
            } else if is_duplicate {
                self.stats.n_duplicate_reads += 1;
                tx_unmapped
                    .send(record.clone())
                    .expect("Failed to send record");
            } else if is_secondary {
                self.stats.n_secondary_reads += 1;
                tx_unmapped
                    .send(record.clone())
                    .expect("Failed to send record");
            } else {
                let ref_seq_id = match record.reference_sequence_id() {
                    Some(Ok(id)) => id,
                    _ => {
                        self.stats.n_unmapped_reads += 1;
                        continue;
                    }
                };

                let is_endogenous = self.endogenous_refseq_ids.contains_key(&ref_seq_id);
                let is_exogenous = self.exogenous_refseq_ids.contains_key(&ref_seq_id);

                if is_endogenous && is_exogenous {
                    self.stats.n_both_genomes += 1;
                    tx_both.send(record.clone()).expect("Failed to send record");
                } else if is_endogenous {
                    self.stats.n_endogenous += 1;
                    tx_endogenous
                        .send(record.clone())
                        .expect("Failed to send record");
                } else if is_exogenous {
                    self.stats.n_exogenous += 1;
                    tx_exogenous
                        .send(record.clone())
                        .expect("Failed to send record");
                }
            }
        }

        // Close the channels
        drop(tx_endogenous);
        drop(tx_exogenous);
        drop(tx_both);
        drop(tx_unmapped);

        // Wait for the writer threads to finish
        endogenous_writer
            .join()
            .expect("Failed to join writer thread");
        exogenous_writer
            .join()
            .expect("Failed to join writer thread");
        both_writer.join().expect("Failed to join writer thread");
        unmapped_writer
            .join()
            .expect("Failed to join writer thread");

        // Print the statistics
        info!("Separation statistics for {}", self.stats.filename);
        info!("{}", self.stats);

        Ok(())
    }

    /// Get the statistics for the split.
    pub fn stats(&self) -> &SplitStats {
        &self.stats
    }

    /// Get the path to the endogenous BAM file.
    pub fn endogenous_file(&self) -> &Path {
        &self.endogenous_bam
    }

    /// Get the path to the exogenous BAM file.
    pub fn exogenous_file(&self) -> &Path {
        &self.exogenous_bam
    }

    /// Get the path to the BAM file containing reads mapping to both genomes.
    pub fn both_file(&self) -> &Path {
        &self.both_bam
    }

    /// Get the path to the BAM file containing unmapped reads.
    pub fn unmapped_file(&self) -> &Path {
        &self.unmapped_bam
    }
}
