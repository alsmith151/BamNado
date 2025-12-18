//! # BAM Modifier Module
//!
//! This module provides functionality for modifying BAM files, such as applying Tn5 shifts
//! to read coordinates and filtering reads based on specific criteria.
//!
//! It is useful for preprocessing BAM files before downstream analysis, ensuring that
//! reads are correctly positioned and filtered.

use crate::bam_utils::get_bam_header;

use crate::read_filter::BamReadFilter;
use anyhow::Context;
use anyhow::Result;

use indicatif::ProgressBar;
use log::warn;

use noodles::bam::r#async::io::{Reader as AsyncReader, Writer as AsyncWriter};
use noodles::bam::bai;
use noodles::core::{Position, Region};

use std::path::PathBuf;
use tokio::fs::File;

/// Modifies BAM files by filtering reads and optionally applying Tn5 shifts.
pub struct BamModifier {
    filepath: PathBuf,
    output: PathBuf,
    filter: BamReadFilter,
    tn5_shift: bool,
}

impl BamModifier {
    /// Creates a new `BamModifier`.
    ///
    /// # Arguments
    ///
    /// * `filepath` - Path to the input BAM file.
    /// * `output` - Path to the output BAM file.
    /// * `filter` - Filter criteria for reads.
    /// * `tn5_shift` - Whether to apply Tn5 shift adjustments.
    pub fn new(filepath: PathBuf, output: PathBuf, filter: BamReadFilter, tn5_shift: bool) -> Self {
        BamModifier {
            filepath,
            output,
            filter,
            tn5_shift,
        }
    }

    /// Runs the modification process asynchronously.
    ///
    /// This reads the input BAM, filters reads, applies shifts if requested,
    /// and writes the result to the output file.
    pub async fn run(&self) -> Result<()> {
        let filepath = self.filepath.clone();
        let outfile = self.output.clone();

        let mut reader = File::open(&filepath).await.map(AsyncReader::new)?;
        let header = reader.read_header().await;

        let header = match header {
            Ok(header) => header,
            Err(_e) => get_bam_header(&filepath)?,
        };

        let index_path = self.filepath.with_extension("bam.bai");
        let index_file = File::open(&index_path).await?;
        let mut index_reader = bai::r#async::io::Reader::new(index_file);
        let index = index_reader.read_index().await?;

        // Make writer
        let mut writer = AsyncWriter::new(File::create(outfile).await?);
        writer.write_header(&header).await?;

        // Get the chromosome sizes
        let chromsizes = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)| (name.to_string(), seq.length().get() as u64))
            .collect::<std::collections::HashMap<_, _>>();

        let query_regions = chromsizes.iter().map(|(name, size)| {
            let start = Position::try_from(1).unwrap();
            let end = Position::try_from(*size as usize).unwrap();
            Region::new(name.to_string(), start..=end)
        });

        let progress = ProgressBar::new(chromsizes.len() as u64);

        for region in query_regions {
            progress.inc(1);
            let mut query = reader.query(&header, &index, &region)?;

            let mut record = noodles::bam::Record::default();
            while query.read_record(&mut record).await? != 0 {
                let is_valid = self.filter.is_valid(&record, Some(&header))?;
                if !is_valid {
                    continue;
                }
                // Apply Tn5 shift if enabled
                if self.tn5_shift {
                    if record.flags().is_properly_segmented() {
                        let shift_values = [4, -5, 5, -4];

                        let reverse = record.flags().is_reverse_complemented();
                        let first_in_template = record.flags().is_first_segment();
                        let tlen = record.template_length();
                        let mut start = record
                            .alignment_start()
                            .context("Missing alignment start")??
                            .get() as i64
                            - 1; // Convert to 0-based

                        start = std::cmp::max(start, 0); // Ensure start is non-negative
                        let mut end = start + record.sequence().len() as i64;

                        let ref_seq_id = record
                            .reference_sequence_id()
                            .context("Missing reference sequence ID")??;
                        let chrom_name = header
                            .reference_sequences()
                            .get_index(ref_seq_id)
                            .context("Missing reference sequence")?
                            .0
                            .to_string();

                        let chromsize = chromsizes
                            .get(&chrom_name)
                            .context("Missing chromosome size")?;

                        let dtlen = match (reverse, first_in_template) {
                            (true, true) => {
                                end += shift_values[1];
                                shift_values[1] - shift_values[0]
                            }
                            (true, false) => {
                                end -= shift_values[2];
                                shift_values[3] - shift_values[2]
                            }
                            (false, true) => {
                                start -= shift_values[3];
                                shift_values[3] - shift_values[2]
                            }
                            (false, false) => {
                                start += shift_values[0];
                                shift_values[1] - shift_values[0]
                            }
                        };

                        // Sanity check coordinates
                        if start >= 0 && end <= *chromsize as i64 {
                            // Update the record - Note: noodles BAM records are immutable
                            // We'll write the record as-is for now

                            let mut aln =
                                noodles::sam::alignment::RecordBuf::try_from_alignment_record(
                                    &header, &record,
                                )?;

                            // Ensure start is within bounds before converting to Position
                            if start >= 0 && start <= *chromsize as i64 {
                                // Position::try_from expects 1-based coordinates, so convert back
                                *aln.alignment_start_mut() =
                                    Some(Position::try_from((start + 1) as usize)?);
                            } else {
                                warn!(
                                    "Adjusted start out of bounds: {start} (chromsize={chromsize})"
                                );
                                // Skip writing this record if start invalid
                                continue;
                            }

                            // Safely adjust template length (tlen is i32)
                            let adjusted_tlen: i32 = if tlen < 0 {
                                // tlen - dtlen
                                tlen.saturating_sub(dtlen as i32)
                            } else {
                                // tlen + dtlen
                                tlen.saturating_add(dtlen as i32)
                            };
                            *aln.template_length_mut() = adjusted_tlen;

                            // Get mate start as signed 0-based coordinate
                            let mate_start = aln
                                .mate_alignment_start()
                                .context("Missing mate alignment start")?
                                .get() as i64
                                - 1; // Convert to 0-based

                            // Compute adjusted mate positions using signed arithmetic
                            // and ensure they are within bounds before converting to usize
                            if !first_in_template && reverse {
                                let adj = mate_start + shift_values[0];
                                if adj >= 0 && adj <= *chromsize as i64 {
                                    // Convert 0-based adj back to 1-based for Position
                                    *aln.mate_alignment_start_mut() =
                                        Some(Position::try_from((adj + 1) as usize)?);
                                } else {
                                    warn!(
                                        "Adjusted mate start out of bounds: {adj} (chromsize={chromsize})"
                                    );
                                }
                            } else if first_in_template && reverse {
                                let adj = mate_start - shift_values[3];
                                if adj >= 0 && adj <= *chromsize as i64 {
                                    // Convert 0-based adj back to 1-based for Position
                                    *aln.mate_alignment_start_mut() =
                                        Some(Position::try_from((adj + 1) as usize)?);
                                } else {
                                    warn!(
                                        "Adjusted mate start out of bounds: {adj} (chromsize={chromsize})"
                                    );
                                }
                            }

                            writer.write_alignment_record(&header, &aln).await?;
                        }
                    } else {
                        warn!("Skipping record with invalid alignment start: {record:?}");
                    }
                } else {
                    // Write the record without TN5 shift
                    writer.write_record(&header, &record).await?;
                }
            }
        }

        progress.finish();
        writer.shutdown().await?;

        Ok(())
    }
}
