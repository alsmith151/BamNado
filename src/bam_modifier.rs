use crate::bam_utils::BamStats;
use crate::bam_utils::get_bam_header;
use crate::genomic_intervals;
use crate::genomic_intervals::{IntervalMaker, Shift};
use crate::read_filter::BamReadFilter;
use anyhow::Context;
use anyhow::Result;
use crossbeam::channel::unbounded;
use indicatif::ProgressBar;
use log::{debug, info, warn};

use noodles::bam::r#async::io::{Reader as AsyncReader, Writer as AsyncWriter};
use noodles::bam::bai;
use noodles::core::{Position, Region};
use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
use noodles::{bam, sam};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::path::PathBuf;

use futures::TryStreamExt;
use tokio::fs::File;

pub struct BamModifier {
    filepath: PathBuf,
    output: PathBuf,
    filter: BamReadFilter,
    tn5_shift: bool,
}

impl BamModifier {
    pub fn new(filepath: PathBuf, output: PathBuf, filter: BamReadFilter, tn5_shift: bool) -> Self {
        BamModifier {
            filepath,
            output,
            filter,
            tn5_shift,
        }
    }

    pub async fn run(&self) -> Result<()> {
        let filepath = self.filepath.clone();
        let outfile = self.output.clone();

        let mut reader = File::open(&filepath).await.map(AsyncReader::new)?;
        let header = reader.read_header().await;

        let header = match header {
            Ok(header) => header,
            Err(e) => get_bam_header(&filepath)?,
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

            // println!("Querying region: {:?}", region);

            let mut query = reader.query(&header, &index, &region)?;
            while let Some(record) = query.try_next().await? {
                let is_valid = self.filter.is_valid(&record, Some(&header))?;
                if !is_valid {
                    continue;
                }

                if self.tn5_shift {
                    let record = record;

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

                            *aln.alignment_start_mut() = Some(Position::try_from(start as usize)?);

                            if tlen < 0 {
                                *aln.template_length_mut() = tlen - dtlen as i32;
                            } else {
                                *aln.template_length_mut() = tlen + dtlen as i32;
                            }

                            let mate_start = aln
                                .mate_alignment_start()
                                .context("Missing mate alignment start")?
                                .get()
                                - 1; // Convert to 0-based

                            if !first_in_template && reverse {
                                *aln.mate_alignment_start_mut() = Some(Position::try_from(
                                    mate_start + shift_values[0] as usize,
                                )?);
                            } else if first_in_template && reverse {
                                *aln.mate_alignment_start_mut() = Some(Position::try_from(
                                    mate_start - shift_values[3] as usize,
                                )?);
                            }


                            writer
                                .write_alignment_record(&header, &aln)
                                .await?;

                        }
                    } else {
                        warn!("Skipping record with invalid alignment start: {:?}", record);
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
