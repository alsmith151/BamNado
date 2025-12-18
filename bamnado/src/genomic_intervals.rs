//! # Genomic Intervals Module
//!
//! This module provides tools for manipulating genomic intervals and coordinates.
//! It includes functionality for:
//! *   Parsing CIGAR strings to determine alignment spans.
//! *   Shifting and truncating genomic coordinates (e.g., for Tn5 adjustments).
//! *   Representing and operating on genomic regions.

use std::cmp::min;
use std::str::FromStr;

use crate::read_filter::BamReadFilter;
use ahash::HashMap;
use anyhow::{Context, Result, anyhow};
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::Cigar;
use noodles::{bam, sam};
use sam::alignment::record::cigar::op::Kind;

/// Returns the alignment span of a CIGAR operation.
/// # Arguments
/// * `cigar` - A reference to a CIGAR operation.
/// # Returns
/// The total length of the alignment span, which is the sum of all match, insertion,
/// and deletion operations in the CIGAR string.
fn alignment_span<C>(cigar: &C) -> Result<usize>
where
    C: Cigar,
{
    cigar.iter().try_fold(0, |acc, op| {
        let op = op.context("Failed to parse CIGAR operation")?;
        match op.kind() {
            Kind::Match
            | Kind::Insertion
            | Kind::SequenceMatch
            | Kind::Pad
            | Kind::SequenceMismatch => Ok(acc + op.len()),
            Kind::Deletion | Kind::Skip => Ok(acc),
            _ => Ok(acc), // Ignore other kinds like SoftClip, HardClip, etc.
        }
    })
}

/// Represents coordinate shifts to be applied to reads.
///
/// The shifts can be applied to the 5' and 3' ends of forward and reverse reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Shift {
    /// Shift for the 5' end of forward reads.
    pub five_prime: i32,
    /// Shift for the 3' end of forward reads.
    pub three_prime: i32,
    /// Shift for the 5' end of reverse reads.
    pub five_prime_reverse: i32,
    /// Shift for the 3' end of reverse reads.
    pub three_prime_reverse: i32,
    index: usize,
}

impl Iterator for Shift {
    type Item = i32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.index {
            0 => {
                self.index += 1;
                Some(self.five_prime)
            }
            1 => {
                self.index += 1;
                Some(self.three_prime)
            }
            2 => {
                self.index += 1;
                Some(self.five_prime_reverse)
            }
            3 => {
                self.index += 1;
                Some(self.three_prime_reverse)
            }
            _ => None,
        }
    }
}

impl FromStr for Shift {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 4 {
            return Err(anyhow!("Shift must be in the format '5,3,5r,3r'"));
        }

        let five_prime = parts[0].parse::<i32>().unwrap_or(0);
        let three_prime = parts[1].parse::<i32>().unwrap_or(0);
        let five_prime_reverse = parts[2].parse::<i32>().unwrap_or(0);
        let three_prime_reverse = parts[3].parse::<i32>().unwrap_or(0);

        Ok(Self {
            five_prime,
            three_prime,
            five_prime_reverse,
            three_prime_reverse,
            index: 0,
        })
    }
}

impl From<[i32; 4]> for Shift {
    fn from(shift: [i32; 4]) -> Self {
        Self {
            five_prime: shift[0],
            three_prime: shift[1],
            five_prime_reverse: shift[2],
            three_prime_reverse: shift[3],
            index: 0,
        }
    }
}
impl From<[usize; 4]> for Shift {
    fn from(shift: [usize; 4]) -> Self {
        Self {
            five_prime: shift[0] as i32,
            three_prime: shift[1] as i32,
            five_prime_reverse: shift[2] as i32,
            three_prime_reverse: shift[3] as i32,
            index: 0,
        }
    }
}

impl From<Vec<i32>> for Shift {
    fn from(shift: Vec<i32>) -> Self {
        if shift.len() != 4 {
            panic!("Shift vector must have exactly 4 elements");
        }
        Self {
            five_prime: shift[0],
            three_prime: shift[1],
            five_prime_reverse: shift[2],
            three_prime_reverse: shift[3],
            index: 0,
        }
    }
}

/// Represents truncation values for reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Truncate {
    five_prime: Option<usize>,
    three_prime: Option<usize>,
}

impl From<(Option<usize>, Option<usize>)> for Truncate {
    fn from(truncate: (Option<usize>, Option<usize>)) -> Self {
        Self {
            five_prime: truncate.0,
            three_prime: truncate.1,
        }
    }
}

impl FromStr for Truncate {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 2 {
            return Err(anyhow!("Truncate must be in the format '5,3'"));
        }

        let five_prime = parts[0].parse::<usize>().ok();
        let three_prime = parts[1].parse::<usize>().ok();

        Ok(Self {
            five_prime,
            three_prime,
        })
    }
}

/// Helper struct for creating genomic intervals from BAM records.
///
/// It handles coordinate calculations, shifting, and truncation.
pub struct IntervalMaker<'a> {
    read: bam::Record,
    header: &'a sam::Header,
    chromsizes: &'a HashMap<usize, u64>,
    filter: &'a BamReadFilter,
    use_fragment: bool,
    shift: Shift,
    truncate: Option<Truncate>,
}

impl<'a> IntervalMaker<'a> {
    /// Creates a new `IntervalMaker`.
    pub fn new(
        read: bam::Record,
        header: &'a sam::Header,
        chromsizes: &'a HashMap<usize, u64>,
        filter: &'a BamReadFilter,
        use_fragment: bool,
        shift: Option<Shift>,
        truncate: Option<Truncate>,
    ) -> Self {
        Self {
            read,
            header,
            chromsizes,
            filter,
            use_fragment,
            shift: shift.map_or_else(|| Shift::from([0, 0, 0, 0]), Shift::from),
            truncate,
        }
    }

    fn coords_fragment(&self) -> Option<(usize, usize, usize)> {
        if !self.read.flags().is_first_segment() {
            return None;
        }

        let template_length = self.read.template_length();

        // If the template length is positive, we use the alignment start of the read.
        // If it's negative, we use the mate alignment start.
        // If it's zero, we return None as this has no valid genomic coordinates.
        if template_length > 0 {
            let start = self
                .read
                .alignment_start()
                .context("Missing alignment start")
                .ok()?
                .ok()?
                .get();
            let end = start + template_length as usize;
            Some((start, end, 0))
        } else if template_length < 0 {
            let start = self
                .read
                .mate_alignment_start()
                .context("Missing mate alignment start")
                .ok()?
                .ok()?
                .get();
            let end = start + template_length.unsigned_abs() as usize;
            Some((start, end, 0))
        } else {
            None
        }
    }

    fn coords_read(&self) -> Option<(usize, usize, usize)> {
        // For single reads, we use the alignment start and calculate the end based on the alignment span.
        // If the read has no alignment start, we return None.
        let start = self
            .read
            .alignment_start()
            .context("Missing alignment start")
            .ok()?
            .ok()?
            .get();
        let alignment_span = alignment_span(&self.read.cigar())
            .context("Failed to calculate alignment span")
            .ok()?;
        let end = start + alignment_span;
        Some((start, end, 0))
    }

    fn coords_shifted(
        &self,
        mut start: usize,
        mut end: usize,
        _dtlen: usize,
    ) -> Result<(usize, usize, i32)> {
        // Changed return type to i32
        let reverse = self.read.flags().is_reverse_complemented();
        let first_in_template = self.read.flags().is_first_segment();

        // Adjust start and end based on the shift values.
        // The logic here depends on whether the read is reverse complemented and whether it is the first
        // segment in the template.
        // The shift values are applied differently based on these conditions.
        // The shift array is expected to have the following meaning:
        // [shift_5_prime, shift_3_prime, shift_5_prime_reverse, shift_3_prime_reverse]
        let dtlen = match (reverse, first_in_template) {
            (true, true) => {
                // Reverse first read: apply shifts correctly for reverse strand
                if self.shift.three_prime >= 0 {
                    start = start.saturating_sub(self.shift.three_prime as usize);
                } else {
                    start = start.saturating_add((-self.shift.three_prime) as usize);
                }
                if self.shift.five_prime_reverse >= 0 {
                    end = end.saturating_add(self.shift.five_prime_reverse as usize);
                } else {
                    end = end.saturating_sub((-self.shift.five_prime_reverse) as usize);
                }
                // For reverse reads, the template length change is the sum of the applied shifts
                self.shift.five_prime_reverse + self.shift.three_prime
            }
            (true, false) => {
                // Reverse second read: use reverse-specific shifts
                if self.shift.three_prime_reverse >= 0 {
                    start = start.saturating_sub(self.shift.three_prime_reverse as usize);
                } else {
                    start = start.saturating_add((-self.shift.three_prime_reverse) as usize);
                }
                if self.shift.five_prime_reverse >= 0 {
                    end = end.saturating_add(self.shift.five_prime_reverse as usize);
                } else {
                    end = end.saturating_sub((-self.shift.five_prime_reverse) as usize);
                }
                self.shift.five_prime_reverse + self.shift.three_prime_reverse
            }
            (false, true) => {
                // Forward first read: direct application of shifts
                if self.shift.five_prime >= 0 {
                    start = start.saturating_sub(self.shift.five_prime as usize);
                } else {
                    start = start.saturating_add((-self.shift.five_prime) as usize);
                }
                if self.shift.three_prime >= 0 {
                    end = end.saturating_add(self.shift.three_prime as usize);
                } else {
                    end = end.saturating_sub((-self.shift.three_prime) as usize);
                }
                self.shift.five_prime + self.shift.three_prime
            }
            (false, false) => {
                // Forward second read: mix of forward and reverse shifts
                if self.shift.five_prime >= 0 {
                    start = start.saturating_sub(self.shift.five_prime as usize);
                } else {
                    start = start.saturating_add((-self.shift.five_prime) as usize);
                }
                if self.shift.three_prime_reverse >= 0 {
                    end = end.saturating_add(self.shift.three_prime_reverse as usize);
                } else {
                    end = end.saturating_sub((-self.shift.three_prime_reverse) as usize);
                }
                self.shift.five_prime + self.shift.three_prime_reverse
            }
        };

        let ref_id = self
            .read
            .reference_sequence_id()
            .context("Missing reference sequence ID")?
            .context("Reference sequence ID is None")?;

        let chromsize = self
            .chromsizes
            .get(&ref_id)
            .ok_or_else(|| anyhow!("No chromosome size found for reference ID {ref_id}"))?;

        // Clamp to valid genomic range
        let start = start.max(1);
        let end = min(*chromsize as usize, end);

        Ok((start, end, dtlen))
    }

    /// Calculates the coordinates of the interval.
    ///
    /// Returns `(start, end, dtlen)` where `dtlen` is the change in template length.
    pub fn coords(&self) -> Option<(usize, usize, i32)> {
        // Changed return type
        match self.filter.is_valid(&self.read, Some(self.header)) {
            Ok(true) => {
                let (start, end, dtlen) = if self.use_fragment {
                    self.coords_fragment()?
                } else {
                    self.coords_read()?
                };

                let (start, end, dtlen) = if self.shift.into_iter().any(|x| x != 0) {
                    self.coords_shifted(start, end, dtlen).ok()?
                } else {
                    (start, end, dtlen as i32)
                };

                let (start, end) = self.truncate(start, end);

                Some((start, end, dtlen))
            }
            _ => None,
        }
    }

    fn truncate(&self, start: usize, end: usize) -> (usize, usize) {
        let mut start = start;
        let mut end = end;

        if self.truncate.is_none() {
            return (start, end);
        }

        // If truncate is set, apply the truncation values.
        // If truncate_5_prime is set, we adjust the start position.
        // If truncate_3_prime is set, we adjust the end position.
        // These values are optional, so we check if they are Some before applying them.
        if let Some(truncate) = self.truncate {
            if let Some(truncate_5_prime) = truncate.five_prime {
                start = start.saturating_add(truncate_5_prime);
            }

            if let Some(truncate_3_prime) = truncate.three_prime {
                end = end.saturating_sub(truncate_3_prime);
            }
        }

        (start, end)
    }

    /// Creates a new BAM record with modified coordinates.
    pub fn record(&self) -> Option<RecordBuf> {
        let (start, _end, dtlen) = self.coords()?;

        let mut new_record = RecordBuf::try_from_alignment_record(self.header, &self.read).ok()?;

        *new_record.alignment_start_mut() = Position::new(start);

        // Update template length based on the delta from coordinate modifications
        let original_template_length = self.read.template_length();

        // Apply the delta to the template length, preserving sign
        let new_template_length = if original_template_length >= 0 {
            original_template_length.saturating_add(dtlen)
        } else {
            original_template_length.saturating_sub(dtlen)
        };

        *new_record.template_length_mut() = new_template_length;

        Some(new_record)
    }
}
