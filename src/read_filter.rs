use ahash::{HashMap, HashSet};
use anyhow::{Context, Result};
use noodles::sam;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::data::field::tag::Tag;
use rust_lapper::Lapper;
use std::fmt::Display;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::bam_utils::CB;

#[derive(Debug)]
pub struct BamReadFilterStats {
    // Total number of reads
    n_total: AtomicU64,
    // Number of reads filtered by proper pair
    n_failed_proper_pair: AtomicU64,
    // Number of reads filtered by mapping quality
    n_failed_mapq: AtomicU64,
    // Number of reads filtered by length
    n_failed_length: AtomicU64,
    // Number of reads filtered by blacklisted locations
    n_failed_blacklist: AtomicU64,
    // Number of reads filtered by barcode
    n_failed_barcode: AtomicU64,
    // Number of reads not in read group
    n_not_in_read_group: AtomicU64,
    // Number of reads filtered by incorrect strand
    n_incorrect_strand: AtomicU64,
}

impl Default for BamReadFilterStats {
    fn default() -> Self {
        Self::new()
    }
}

impl BamReadFilterStats {
    pub fn new() -> Self {
        Self {
            n_total: AtomicU64::new(0),
            n_failed_proper_pair: AtomicU64::new(0),
            n_failed_mapq: AtomicU64::new(0),
            n_failed_length: AtomicU64::new(0),
            n_failed_blacklist: AtomicU64::new(0),
            n_failed_barcode: AtomicU64::new(0),
            n_not_in_read_group: AtomicU64::new(0),
            n_incorrect_strand: AtomicU64::new(0),
        }
    }

    pub fn snapshot(&self) -> BamReadFilterStatsSnapshot {
        BamReadFilterStatsSnapshot {
            n_total: self.n_total.load(Ordering::Relaxed),
            n_failed_proper_pair: self.n_failed_proper_pair.load(Ordering::Relaxed),
            n_failed_mapq: self.n_failed_mapq.load(Ordering::Relaxed),
            n_failed_length: self.n_failed_length.load(Ordering::Relaxed),
            n_failed_blacklist: self.n_failed_blacklist.load(Ordering::Relaxed),
            n_failed_barcode: self.n_failed_barcode.load(Ordering::Relaxed),
            n_not_in_read_group: self.n_not_in_read_group.load(Ordering::Relaxed),
            n_incorrect_strand: self.n_incorrect_strand.load(Ordering::Relaxed),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BamReadFilterStatsSnapshot {
    n_total: u64,
    n_failed_proper_pair: u64,
    n_failed_mapq: u64,
    n_failed_length: u64,
    n_failed_blacklist: u64,
    n_failed_barcode: u64,
    n_not_in_read_group: u64,
    n_incorrect_strand: u64,
}

impl Display for BamReadFilterStatsSnapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f)?;
        writeln!(f, "Total reads: {}", self.n_total)?;
        writeln!(f, "Failed proper pair: {}", self.n_failed_proper_pair)?;
        writeln!(f, "Failed mapping quality: {}", self.n_failed_mapq)?;
        writeln!(f, "Failed length: {}", self.n_failed_length)?;
        writeln!(f, "Failed blacklist: {}", self.n_failed_blacklist)?;
        writeln!(f, "Failed barcode: {}", self.n_failed_barcode)?;
        writeln!(f, "Not in read group: {}", self.n_not_in_read_group)?;
        writeln!(f, "Incorrect strand: {}", self.n_incorrect_strand)?;
        writeln!(
            f,
            "Filtered reads: {}",
            self.n_failed_proper_pair
                + self.n_failed_mapq
                + self.n_failed_length
                + self.n_failed_blacklist
                + self.n_failed_barcode
                + self.n_not_in_read_group
                + self.n_incorrect_strand
        )?;
        Ok(())
    }
}

impl BamReadFilterStatsSnapshot {
    pub fn n_reads_after_filtering(&self) -> u64 {
        self.n_total
            - (self.n_failed_proper_pair
                + self.n_failed_mapq
                + self.n_failed_length
                + self.n_failed_blacklist
                + self.n_failed_barcode
                + self.n_not_in_read_group
                + self.n_incorrect_strand)
    }
}

/// A filter for BAM reads.
/// Set the minimum mapping quality, minimum and maximum read length, blacklisted locations, and whitelisted barcodes.
/// The filter is applied to each read in the BAM file.
#[derive(Debug, Clone)]
pub struct BamReadFilter {
    // Whether to select only a specific strand of reads
    strand: bio_types::strand::Strand,
    // Properly paired reads only
    proper_pair: bool,
    // Minimum mapping quality
    min_mapq: u8,
    // Minimum read length
    min_length: u32,
    // Maximum read length
    max_length: u32,
    // Blacklisted locations (chromosome (reference number) -> Lapper)
    blacklisted_locations: Option<HashMap<usize, Lapper<usize, u32>>>,
    // Whitelisted barcodes (cell barcodes to keep)
    whitelisted_barcodes: Option<HashSet<String>>,
    // Read group to keep
    read_group: Option<String>,
    // Statistics for the filtering process
    stats: Arc<BamReadFilterStats>,
}

impl Display for BamReadFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "\tOnly allowing reads with strand: {}", self.strand)?;
        writeln!(f, "\tProper pair: {}", self.proper_pair)?;
        writeln!(f, "\tMinimum mapping quality: {}", self.min_mapq)?;
        writeln!(f, "\tMinimum read length: {}", self.min_length)?;
        writeln!(f, "\tMaximum read length: {}", self.max_length)?;

        match &self.blacklisted_locations {
            Some(blacklisted_locations) => {
                writeln!(
                    f,
                    "\tNumber of Blacklisted locations: {}",
                    blacklisted_locations
                        .values()
                        .map(|v| v.len())
                        .sum::<usize>()
                )?;
            }
            None => {
                writeln!(f, "\tNumber of Blacklisted locations: 0")?;
            }
        }

        match &self.whitelisted_barcodes {
            Some(whitelisted_barcodes) => {
                writeln!(
                    f,
                    "\tNumber of Whitelisted barcodes: {}",
                    whitelisted_barcodes.len()
                )?;
            }
            None => {
                writeln!(f, "\tNumber of Whitelisted barcodes: 0")?;
            }
        }
        Ok(())
    }
}

impl Default for BamReadFilter {
    fn default() -> Self {
        Self::new(
            bio_types::strand::Strand::Unknown,
            true,
            Some(0),
            Some(0),
            Some(1000),
            None,
            None,
            None,
        )
    }
}

impl BamReadFilter {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        strand: bio_types::strand::Strand,
        proper_pair: bool,
        min_mapq: Option<u8>,
        min_length: Option<u32>,
        max_length: Option<u32>,
        read_group: Option<String>,
        blacklisted_locations: Option<HashMap<usize, Lapper<usize, u32>>>,
        whitelisted_barcodes: Option<HashSet<String>>,
    ) -> Self {
        let min_mapq = min_mapq.unwrap_or(0);
        let min_length = min_length.unwrap_or(0);
        let max_length = max_length.unwrap_or(u32::MAX);

        Self {
            strand,
            proper_pair,
            min_mapq,
            min_length,
            max_length,
            blacklisted_locations,
            whitelisted_barcodes,
            read_group,
            stats: Arc::new(BamReadFilterStats::new()),
        }
    }

    pub fn is_valid<R>(&self, alignment: &R, header: Option<&sam::Header>) -> Result<bool>
    where
        R: sam::alignment::Record,
    {
        // Update total count.
        self.stats.n_total.fetch_add(1, Ordering::Relaxed);

        let flags = match alignment.flags() {
            Ok(flags) => flags,
            Err(_) => {
                self.stats
                    .n_failed_proper_pair
                    .fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }
        };

        // Filter by strand.
        match flags.is_reverse_complemented() {
            true => {
                // Read is reverse stranded. If the filter is set for only forward reads we need to return false.
                if self.strand == bio_types::strand::Strand::Forward {
                    self.stats
                        .n_incorrect_strand
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
            false => {
                // Read is forward stranded. If the filter is set for only reverse reads we need to return false.
                if self.strand == bio_types::strand::Strand::Reverse {
                    self.stats
                        .n_incorrect_strand
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
        }

        // Filter by proper pair.
        if self.proper_pair && !flags.is_properly_segmented() {
            return Ok(false);
        }

        // Filter by unmapped reads.
        if flags.is_unmapped() {
            self.stats.n_failed_mapq.fetch_add(1, Ordering::Relaxed);
            return Ok(false);
        }

        // Filter by mapping quality.
        match alignment.mapping_quality() {
            Some(Ok(mapping_quality)) => {
                if mapping_quality.get() < self.min_mapq {
                    self.stats.n_failed_mapq.fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
            None => {
                // STAR mapping quality does not match the standard SAM format will just assume it is maximum.
                // This is a workaround for STAR mappings that do not provide mapping quality.
            }
            _ => {
                self.stats.n_failed_mapq.fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }
        }

        // Filter by read length.
        let alignment_length = alignment.sequence().len();
        if alignment_length < self.min_length as usize
            || alignment_length > self.max_length as usize
        {
            self.stats.n_failed_length.fetch_add(1, Ordering::Relaxed);
            return Ok(false);
        }

        let header = header.expect("No header provided");
        let chrom_id = alignment
            .reference_sequence_id(header)
            .context("Failed to get reference sequence ID")??;

        let start = match alignment.alignment_start() {
            Some(Ok(start)) => start.get(),
            _ => {
                self.stats.n_failed_mapq.fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }
        };
        let end = start + alignment_length;

        // Filter by blacklisted locations.
        if let Some(blacklisted_locations) = &self.blacklisted_locations {
            if let Some(blacklist) = blacklisted_locations.get(&chrom_id) {
                if blacklist.count(start, end) > 0 {
                    self.stats
                        .n_failed_blacklist
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
        }

        // Filter by whitelisted barcodes.
        if let Some(barcodes) = &self.whitelisted_barcodes {
            let barcode = get_cell_barcode(alignment);
            match barcode {
                Some(barcode) => {
                    if !barcodes.contains(&barcode) {
                        self.stats.n_failed_barcode.fetch_add(1, Ordering::Relaxed);
                        return Ok(false);
                    }
                }
                None => {
                    self.stats.n_failed_barcode.fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
        }

        // Filter by read group.
        if let Some(read_group) = &self.read_group {
            let rg = noodles::sam::alignment::record::data::field::tag::Tag::READ_GROUP;
            let data = alignment.data();
            let read_group_value = data.get(&rg).context("Failed to get read group")??;

            match read_group_value {
                Value::String(value) => {
                    if value != read_group {
                        self.stats
                            .n_not_in_read_group
                            .fetch_add(1, Ordering::Relaxed);
                        return Ok(false);
                    }
                }
                _ => {
                    self.stats
                        .n_not_in_read_group
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
        }

        Ok(true)
    }

    pub fn stats(&self) -> BamReadFilterStatsSnapshot {
        self.stats.snapshot()
    }
}

/// Get the cell barcode from a BAM alignment.
fn get_cell_barcode<R>(alignment: &R) -> Option<String>
where
    R: sam::alignment::Record,
{
    let tags = alignment.data();
    let cell_barcode_tag = Tag::from(CB);
    let tag_value = tags.get(&cell_barcode_tag);
    if let Some(Ok(Value::String(barcode))) = tag_value {
        Some(barcode.to_string())
    } else {
        None
    }
}
