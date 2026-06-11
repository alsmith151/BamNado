//! Read-level manipulation operations.
//!
//! This module is a home for operations that transform individual reads or extract
//! sub-sequences from them. `SoftclipExtractor` is the first; others can be added here.

use crate::read_filter::BamReadFilter;
use anyhow::{Context, Result};
use log::info;
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::{bam, fastq};
use noodles::fastq::record::Definition;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

/// Which soft-clipped ends to extract from each read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ClipEnd {
    #[value(name = "both")]
    Both,
    #[value(name = "5p")]
    FivePrime,
    #[value(name = "3p")]
    ThreePrime,
}

/// Extracts soft-clipped bases from aligned BAM reads and writes them as FASTQ.
///
/// Designed for workflows where reads are aligned to a primary genome and the
/// unaligned (soft-clipped) portions need re-alignment to a secondary genome
/// (e.g. spike-in separation, multi-genome analysis).
pub struct SoftclipExtractor {
    bam_path: PathBuf,
    output_path: PathBuf,
    min_clip_len: usize,
    clip_end: ClipEnd,
    filter: BamReadFilter,
}

impl SoftclipExtractor {
    pub fn new(
        bam_path: PathBuf,
        output_path: PathBuf,
        min_clip_len: usize,
        clip_end: ClipEnd,
        filter: BamReadFilter,
    ) -> Self {
        Self {
            bam_path,
            output_path,
            min_clip_len,
            clip_end,
            filter,
        }
    }

    /// Run extraction. Returns (reads_seen, clips_written).
    pub fn run(&self) -> Result<(u64, u64)> {
        let mut reader = bam::io::reader::Builder
            .build_from_path(&self.bam_path)
            .context("Failed to open BAM file")?;
        let header = reader.read_header().context("Failed to read BAM header")?;

        let fq_file =
            BufWriter::new(File::create(&self.output_path).context("Failed to create output FASTQ")?);
        let mut writer = fastq::io::Writer::new(fq_file);

        let (mut reads_seen, mut clips_written) = (0u64, 0u64);

        for result in reader.records() {
            let record = result.context("Failed to read BAM record")?;
            reads_seen += 1;

            if !self.filter.is_valid(&record, Some(&header))? {
                continue;
            }

            let seq_len = record.sequence().len();
            if seq_len == 0 {
                continue;
            }

            let read_name = record
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_default();

            let (five_clip, three_clip) =
                find_soft_clips(&record.cigar()).context("Failed to parse CIGAR")?;

            for (range_opt, suffix) in [(five_clip, "_5p"), (three_clip, "_3p")] {
                let emit = match self.clip_end {
                    ClipEnd::Both => range_opt.is_some(),
                    ClipEnd::FivePrime => suffix == "_5p" && range_opt.is_some(),
                    ClipEnd::ThreePrime => suffix == "_3p" && range_opt.is_some(),
                };
                if !emit {
                    continue;
                }
                let (start, end) = range_opt.unwrap();
                if end - start < self.min_clip_len {
                    continue;
                }

                let (seq, qual) = extract_clip(&record, start, end);
                let def = Definition::new(format!("{}{}", read_name, suffix), "");
                let fq_record = fastq::Record::new(def, seq, qual);
                writer
                    .write_record(&fq_record)
                    .context("Failed to write FASTQ record")?;
                clips_written += 1;
            }
        }

        info!("Reads processed: {reads_seen}, clips written: {clips_written}");
        Ok((reads_seen, clips_written))
    }
}

/// Walk the CIGAR tracking query position to find soft-clip ranges.
///
/// Uses cumulative query-position tracking so that insertions in the aligned
/// region correctly shift the trailing clip start — rather than the fragile
/// `seq_len - trailing_clip_len` shortcut that breaks on non-standard CIGARs.
fn find_soft_clips<C>(cigar: &C) -> Result<(Option<(usize, usize)>, Option<(usize, usize)>)>
where
    C: Cigar,
{
    let mut query_pos: usize = 0;
    let mut seen_aligned = false;
    let mut five_clip: Option<(usize, usize)> = None;
    let mut three_clip: Option<(usize, usize)> = None;

    for result in cigar.iter() {
        let op = result.context("Failed to parse CIGAR operation")?;
        let len = op.len();

        match op.kind() {
            Kind::SoftClip => {
                let range = (query_pos, query_pos + len);
                if !seen_aligned {
                    five_clip = Some(range);
                } else {
                    three_clip = Some(range);
                }
                query_pos += len;
            }
            Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch => {
                seen_aligned = true;
                query_pos += len;
            }
            // HardClip: at CIGAR boundary but NOT in query sequence — do not advance query_pos.
            // Deletion/Skip/Pad: consume reference only — do not advance query_pos.
            Kind::HardClip | Kind::Deletion | Kind::Skip | Kind::Pad => {}
        }
    }

    Ok((five_clip, three_clip))
}

/// Extract the slice [start, end) of sequence and quality from a BAM record.
///
/// Quality scores are converted from raw Phred (0-93) to FASTQ ASCII (+33).
/// If quality is absent (stored as 0xFF in BAM), fills with b'!' (Phred 0).
fn extract_clip(record: &bam::Record, start: usize, end: usize) -> (Vec<u8>, Vec<u8>) {
    let len = end - start;

    let seq: Vec<u8> = record.sequence().iter().skip(start).take(len).collect();

    let qual_raw: Vec<u8> = record
        .quality_scores()
        .iter()
        .skip(start)
        .take(len)
        .collect();

    let qual = if qual_raw.is_empty() || qual_raw.first() == Some(&0xFF) {
        vec![b'!'; len]
    } else {
        qual_raw.iter().map(|&q| q.saturating_add(33)).collect()
    };

    (seq, qual)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record_buf::Cigar as CigarBuf;
    use tempfile::TempDir;

    fn test_bam() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("Failed to get workspace root")
            .join("test/data/test.bam")
    }

    fn make_cigar(ops: &[(Kind, usize)]) -> CigarBuf {
        CigarBuf::from(ops.iter().map(|&(k, l)| Op::new(k, l)).collect::<Vec<_>>())
    }

    #[test]
    fn test_find_clips_both_ends() {
        let cigar = make_cigar(&[(Kind::SoftClip, 10), (Kind::Match, 80), (Kind::SoftClip, 10)]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, Some((0, 10)));
        assert_eq!(three, Some((90, 100)));
    }

    #[test]
    fn test_find_clips_five_prime_only() {
        let cigar = make_cigar(&[(Kind::SoftClip, 5), (Kind::Match, 95)]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, Some((0, 5)));
        assert_eq!(three, None);
    }

    #[test]
    fn test_find_clips_three_prime_only() {
        let cigar = make_cigar(&[(Kind::Match, 95), (Kind::SoftClip, 5)]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, None);
        assert_eq!(three, Some((95, 100)));
    }

    #[test]
    fn test_find_clips_no_clips() {
        let cigar = make_cigar(&[(Kind::Match, 100)]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, None);
        assert_eq!(three, None);
    }

    #[test]
    fn test_find_clips_insertion_shifts_trailing() {
        // 10S 5I 75M 10S — total query length = 10+5+75+10 = 100
        // insertion (5I) advances query_pos so trailing clip starts at 10+5+75 = 90
        let cigar = make_cigar(&[
            (Kind::SoftClip, 10),
            (Kind::Insertion, 5),
            (Kind::Match, 75),
            (Kind::SoftClip, 10),
        ]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, Some((0, 10)));
        assert_eq!(three, Some((90, 100)));
    }

    #[test]
    fn test_find_clips_deletion_does_not_shift() {
        // 10S 5D 85M 10S — deletion does NOT consume query bases
        // total query length = 10+85+10 = 105 (deletion adds 5 to reference, not query)
        let cigar = make_cigar(&[
            (Kind::SoftClip, 10),
            (Kind::Deletion, 5),
            (Kind::Match, 85),
            (Kind::SoftClip, 10),
        ]);
        let (five, three) = find_soft_clips(&cigar).unwrap();
        assert_eq!(five, Some((0, 10)));
        assert_eq!(three, Some((95, 105)));
    }

    #[test]
    fn test_extractor_produces_output() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output = temp_dir.path().join("clips.fastq");

        let extractor = SoftclipExtractor::new(
            test_bam(),
            output.clone(),
            1,
            ClipEnd::Both,
            BamReadFilter::default(),
        );

        let (reads_seen, _) = extractor.run().expect("Extractor failed");
        assert!(reads_seen > 0);
        assert!(output.exists());
    }

    #[test]
    fn test_min_clip_len_filter() {
        let temp_dir = TempDir::new().expect("Failed to create temp dir");
        let output_low = temp_dir.path().join("clips_low.fastq");
        let output_high = temp_dir.path().join("clips_high.fastq");

        let ext_low = SoftclipExtractor::new(
            test_bam(),
            output_low.clone(),
            1,
            ClipEnd::Both,
            BamReadFilter::default(),
        );
        let ext_high = SoftclipExtractor::new(
            test_bam(),
            output_high.clone(),
            10_000,
            ClipEnd::Both,
            BamReadFilter::default(),
        );

        let (_, clips_low) = ext_low.run().unwrap();
        let (_, clips_high) = ext_high.run().unwrap();
        assert!(clips_low >= clips_high);
    }
}
