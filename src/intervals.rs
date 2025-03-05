use std::cmp::min;

use ahash::HashMap;
use anyhow::{Context, Result};
use noodles::core::Position;

use noodles::{bam, sam};
use noodles::sam::alignment::RecordBuf;
use crate::filter::BamReadFilter;

pub struct IntervalMaker<'a> {
    // The read to create the interval from
    read: bam::Record,

    // The header of the BAM file
    header: &'a sam::Header,

    // The chromosome dictionary (chromosome -> length)
    chromsizes: &'a HashMap<usize, u64>,

    // Interval Filterer
    filter: &'a BamReadFilter,

    // Use the fragment to create the interval instead of the read
    use_fragment: bool,

    // Shift the interval by this amount
    shift: [usize; 4],
}



impl IntervalMaker<'_> {
    pub fn new<'a>(
        read: bam::Record,
        header: &'a sam::Header,
        chromsizes: &'a HashMap<usize, u64>,
        filter: &'a BamReadFilter,
        use_fragment: bool,
        shift: Option<[usize; 4]>,
    ) -> IntervalMaker<'a> {
        let shift = match shift {
            Some(shift) => shift,
            None => [0, 0, 0, 0],
        };

        IntervalMaker {
            read,
            header,
            chromsizes,
            filter,
            use_fragment,
            shift,
        }
    }

    fn coords_fragment(&self) -> Option<(usize, usize, usize)> {
        // Only compute coordinates for the first segment to avoid duplicating paired-end fragments.
        if !self.read.flags().is_first_segment() {
            return None;
        }
    
        let template_length = self.read.template_length();
        if template_length > 0 {
            // For a positive template length, the first read is leftmost.
            // Use its alignment start as the fragment start.

            let start = self.read
                .alignment_start()
                .expect("No start position")
                .expect("Error with getting start")
                .get();
            let end = start + template_length as usize;
            Some((start, end, 0))
        } else if template_length < 0 {
            // For a negative template length, the mate is leftmost.
            let start = self.read
                .mate_alignment_start()
                .expect("No start position")
                .expect("Error with getting start")
                .get();
            let end = start + template_length.abs() as usize;
            Some((start, end, 0))
        } else {
            // Template length of zero indicates an invalid or degenerate fragment.
            None
        }
    }

    fn coords_read(&self) -> Option<(usize, usize, usize)> {
            let read = &self.read;
            let start = read
                .alignment_start()
                .expect("Failed to get start")
                .expect("Error with getting start")
                .get();
            let end = start + read.sequence().len();
            Some((start, end, 0))
    }

    /// Get the shifted coordinates of the read
    /// Returns (start, end, dtlen)
    /// If the read is filtered, return None
    fn coords_shifted(&self) -> Option<(usize, usize, usize)> {

        let reverse = self.read.flags().is_reverse_complemented();
        let first_in_template = self.read.flags().is_first_segment();

        let mut start = self
            .read
            .alignment_start()
            .expect("Failed to get start")
            .expect("Error with getting start")
            .get();
        let mut end = start as usize + self.read.sequence().len();

        let dtlen = match (reverse, first_in_template) {
            (true, true) => {
                end += self.shift[1];
                self.shift[1] - self.shift[0]
            }
            (true, false) => {
                end -= self.shift[2];
                self.shift[3] - self.shift[2]
            }
            (false, true) => {
                start -= self.shift[3];
                self.shift[3] - self.shift[2]
            }
            (false, false) => {
                start += self.shift[0];
                self.shift[1] - self.shift[0]
            }
        };

        let chromsize = self.chromsizes.get(
            &self
                .read
                .reference_sequence_id()
                .expect("No reference sequence ID")
                .expect("Failed to get reference sequence ID"),
        );

        match chromsize {
            Some(chromsize) => {
                let start = min(1, start);
                let end = min(*chromsize as usize, end);
                Some((start, end, dtlen))
            }
            None => None,
        }
    }

    pub fn coords(&self) -> Option<(usize, usize)> {
        match self.filter.is_valid(&self.read, Some(self.header)) {
            Ok(true) => {
                let coord = if self.use_fragment {
                    self.coords_fragment()
                } else if self.shift.iter().any(|&x| x != 0) {
                    self.coords_shifted()
                } else {
                    self.coords_read()
                };

                match coord {
                    None => {
                        None
                    }
                    Some((start, end, _)) => {
                        Some((start, end))
                    }
                }
            }
            _ => None,
        }
    }

    pub fn record(&self) -> Option<RecordBuf> {
        let need_to_shift = self.shift.iter().any(|&x| x != 0);

        let coordinates = if self.use_fragment {
            self.coords_fragment()
        } else if need_to_shift {
            self.coords_shifted()
        } else {
            self.coords_read()
        };

        match coordinates {
            None => None,
            Some((start, end, dtlen)) => {
                let mut record_new = RecordBuf::try_from_alignment_record(self.header, &self.read)
                    .expect("Failed to convert record");

                *record_new.alignment_start_mut() = Position::new(start);
                *record_new.template_length_mut() = match self.read.template_length() > 0 {
                    true => self.read.template_length() + dtlen as i32,
                    false => self.read.template_length() - dtlen as i32,
                };

                Some(record_new)
            }
        }
    }

}
