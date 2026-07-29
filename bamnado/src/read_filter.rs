//! # Read Filter Module
//!
//! This module defines structures and logic for filtering BAM reads. It allows users to
//! include or exclude reads based on:
//! *   Mapping quality (MAPQ).
//! *   SAM flags (e.g., proper pair, unmapped, duplicate).
//! *   Read length.
//! *   Genomic blacklists.
//! *   Cell barcodes and read groups.
//! *   Strand orientation.

use ahash::{HashMap, HashSet};
use anyhow::{Context, Result};
use comfy_table::{
    Cell, ContentArrangement, Table, modifiers::UTF8_ROUND_CORNERS, presets::UTF8_FULL,
};
use noodles::sam;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::data::field::tag::Tag;
use rust_lapper::Lapper;
use std::fmt::Display;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::bam_utils::CB;

/// Statistics for the read filtering process.
///
/// Tracks the number of reads filtered by each criterion.
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
    // Number of reads filtered by tag
    n_failed_tag_filter: AtomicU64,
    // Number of reads filtered by fragment length (insert size / TLEN)
    n_failed_fragment_length: AtomicU64,
    // Number of reads filtered as PCR/optical duplicates
    n_failed_duplicate: AtomicU64,
    // Number of reads filtered as secondary alignments
    n_failed_secondary: AtomicU64,
    // Number of reads filtered as supplementary alignments
    n_failed_supplementary: AtomicU64,
}

impl Default for BamReadFilterStats {
    fn default() -> Self {
        Self::new()
    }
}

impl BamReadFilterStats {
    /// Creates a new `BamReadFilterStats` with all counters set to zero.
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
            n_failed_tag_filter: AtomicU64::new(0),
            n_failed_fragment_length: AtomicU64::new(0),
            n_failed_duplicate: AtomicU64::new(0),
            n_failed_secondary: AtomicU64::new(0),
            n_failed_supplementary: AtomicU64::new(0),
        }
    }

    /// Takes a snapshot of the current statistics.
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
            n_failed_tag_filter: self.n_failed_tag_filter.load(Ordering::Relaxed),
            n_failed_fragment_length: self.n_failed_fragment_length.load(Ordering::Relaxed),
            n_failed_duplicate: self.n_failed_duplicate.load(Ordering::Relaxed),
            n_failed_secondary: self.n_failed_secondary.load(Ordering::Relaxed),
            n_failed_supplementary: self.n_failed_supplementary.load(Ordering::Relaxed),
        }
    }
}

/// A snapshot of the read filter statistics at a specific point in time.
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
    n_failed_tag_filter: u64,
    n_failed_fragment_length: u64,
    n_failed_duplicate: u64,
    n_failed_secondary: u64,
    n_failed_supplementary: u64,
}

impl Display for BamReadFilterStatsSnapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let rows = self.stage_rows();
        let mut table = Table::new();
        table
            .load_preset(UTF8_FULL)
            .apply_modifier(UTF8_ROUND_CORNERS)
            .set_content_arrangement(ContentArrangement::Dynamic)
            .set_header(vec!["stage", "remain", "drop", "% total"]);

        for (label, remain, dropped, pct) in rows {
            let dropped = if dropped == 0 {
                "-".to_string()
            } else {
                format!("-{}", Self::format_count(dropped))
            };

            table.add_row(vec![
                Cell::new(label),
                Cell::new(Self::format_count(remain))
                    .set_alignment(comfy_table::CellAlignment::Right),
                Cell::new(dropped).set_alignment(comfy_table::CellAlignment::Right),
                Cell::new(format!("{pct:.1}%")).set_alignment(comfy_table::CellAlignment::Right),
            ]);
        }

        writeln!(f, "Read Filtering Funnel")?;
        write!(f, "{table}")?;

        Ok(())
    }
}

impl BamReadFilterStatsSnapshot {
    fn format_count(count: u64) -> String {
        let digits = count.to_string();
        let mut out = String::with_capacity(digits.len() + digits.len() / 3);
        for (i, ch) in digits.chars().rev().enumerate() {
            if i > 0 && i % 3 == 0 {
                out.push(',');
            }
            out.push(ch);
        }
        out.chars().rev().collect()
    }

    /// Returns the number of reads filtered out.
    pub fn n_filtered(&self) -> u64 {
        self.n_total - self.n_reads_after_filtering()
    }

    /// Returns the percentage of reads that passed filtering.
    pub fn pass_rate(&self) -> f64 {
        if self.n_total == 0 {
            0.0
        } else {
            (self.n_reads_after_filtering() as f64 / self.n_total as f64) * 100.0
        }
    }

    fn pct_of_total(&self, count: u64) -> f64 {
        if self.n_total == 0 {
            0.0
        } else {
            (count as f64 / self.n_total as f64) * 100.0
        }
    }

    /// Returns rows for a staged read-filtering summary.
    pub fn stage_rows(&self) -> Vec<(&'static str, u64, u64, f64)> {
        let mut remaining = self.n_total;
        let mut rows = vec![("Initial reads", remaining, 0, 100.0)];

        let failures = [
            ("Not a duplicate", self.n_failed_duplicate),
            ("Not secondary", self.n_failed_secondary),
            ("Not supplementary", self.n_failed_supplementary),
            ("Pass strand filter", self.n_incorrect_strand),
            ("Pass proper-pair filter", self.n_failed_proper_pair),
            ("Pass minimum MAPQ", self.n_failed_mapq),
            ("Pass read-length filter", self.n_failed_length),
            ("Pass fragment-length filter", self.n_failed_fragment_length),
            ("Outside blacklist", self.n_failed_blacklist),
            ("Pass barcode allowlist", self.n_failed_barcode),
            ("Pass read-group filter", self.n_not_in_read_group),
            ("Pass tag filter", self.n_failed_tag_filter),
        ];

        for (label, dropped) in failures {
            if dropped > 0 {
                remaining -= dropped;
                rows.push((label, remaining, dropped, self.pct_of_total(remaining)));
            }
        }

        rows.push((
            "Reads retained",
            self.n_reads_after_filtering(),
            0,
            self.pass_rate(),
        ));
        rows
    }

    /// Returns the number of reads remaining after filtering.
    pub fn n_reads_after_filtering(&self) -> u64 {
        self.n_total
            - (self.n_failed_proper_pair
                + self.n_failed_mapq
                + self.n_failed_length
                + self.n_failed_blacklist
                + self.n_failed_barcode
                + self.n_not_in_read_group
                + self.n_incorrect_strand
                + self.n_failed_tag_filter
                + self.n_failed_fragment_length
                + self.n_failed_duplicate
                + self.n_failed_secondary
                + self.n_failed_supplementary)
    }

    /// Returns the total number of reads processed.
    pub fn n_total(&self) -> u64 {
        self.n_total
    }

    /// Returns the number of reads filtered by mapping quality.
    pub fn n_failed_mapq(&self) -> u64 {
        self.n_failed_mapq
    }

    /// Returns the number of reads filtered by length.
    pub fn n_failed_length(&self) -> u64 {
        self.n_failed_length
    }

    /// Returns the number of reads filtered by incorrect strand.
    pub fn n_incorrect_strand(&self) -> u64 {
        self.n_incorrect_strand
    }

    /// Returns the number of reads filtered by proper pair criterion.
    pub fn n_failed_proper_pair(&self) -> u64 {
        self.n_failed_proper_pair
    }

    /// Returns the number of reads filtered by blacklist.
    pub fn n_failed_blacklist(&self) -> u64 {
        self.n_failed_blacklist
    }

    /// Returns the number of reads filtered by barcode.
    pub fn n_failed_barcode(&self) -> u64 {
        self.n_failed_barcode
    }

    /// Returns the number of reads not in the specified read group.
    pub fn n_not_in_read_group(&self) -> u64 {
        self.n_not_in_read_group
    }

    /// Returns the number of reads filtered by tag filter.
    pub fn n_failed_tag_filter(&self) -> u64 {
        self.n_failed_tag_filter
    }

    /// Returns the number of reads filtered by fragment length.
    pub fn n_failed_fragment_length(&self) -> u64 {
        self.n_failed_fragment_length
    }

    /// Returns the number of reads filtered as PCR/optical duplicates.
    pub fn n_failed_duplicate(&self) -> u64 {
        self.n_failed_duplicate
    }

    /// Returns the number of reads filtered as secondary alignments.
    pub fn n_failed_secondary(&self) -> u64 {
        self.n_failed_secondary
    }

    /// Returns the number of reads filtered as supplementary alignments.
    pub fn n_failed_supplementary(&self) -> u64 {
        self.n_failed_supplementary
    }
}

/// A filter for BAM reads.
///
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
    // SAM tag to filter by
    filter_tag: Option<String>,
    // Value for the filter tag
    filter_tag_value: Option<String>,
    // Minimum fragment length (template length / insert size)
    min_fragment_length: Option<u32>,
    // Maximum fragment length (template length / insert size)
    max_fragment_length: Option<u32>,
    // Exclude PCR/optical duplicate reads
    exclude_duplicates: bool,
    // Exclude secondary alignments
    exclude_secondary: bool,
    // Exclude supplementary alignments
    exclude_supplementary: bool,
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
        writeln!(f, "\tExclude duplicates: {}", self.exclude_duplicates)?;
        writeln!(f, "\tExclude secondary: {}", self.exclude_secondary)?;
        writeln!(f, "\tExclude supplementary: {}", self.exclude_supplementary)?;

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

        match (&self.filter_tag, &self.filter_tag_value) {
            (Some(tag), Some(value)) => {
                writeln!(f, "\tFilter tag: {} = {}", tag, value)?;
            }
            _ => {
                writeln!(f, "\tFilter tag: None")?;
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
            None,
            None,
            None,
            None,
        )
    }
}

impl BamReadFilter {
    /// Creates a new `BamReadFilter`.
    ///
    /// # Arguments
    ///
    /// * `strand` - The strand to keep.
    /// * `proper_pair` - Whether to keep only properly paired reads.
    /// * `min_mapq` - Minimum mapping quality.
    /// * `min_length` - Minimum read length.
    /// * `max_length` - Maximum read length.
    /// * `read_group` - Read group to keep.
    /// * `blacklisted_locations` - Genomic regions to exclude.
    /// * `whitelisted_barcodes` - Cell barcodes to keep.
    /// * `filter_tag` - SAM tag to filter by.
    /// * `filter_tag_value` - Value for the filter tag.
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
        filter_tag: Option<String>,
        filter_tag_value: Option<String>,
        min_fragment_length: Option<u32>,
        max_fragment_length: Option<u32>,
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
            filter_tag,
            filter_tag_value,
            min_fragment_length,
            max_fragment_length,
            exclude_duplicates: false,
            exclude_secondary: false,
            exclude_supplementary: false,
            stats: Arc::new(BamReadFilterStats::new()),
        }
    }

    /// Exclude PCR/optical duplicate reads (SAM flag 0x400).
    pub fn with_exclude_duplicates(mut self, v: bool) -> Self {
        self.exclude_duplicates = v;
        self
    }

    /// Exclude secondary alignments (SAM flag 0x100).
    pub fn with_exclude_secondary(mut self, v: bool) -> Self {
        self.exclude_secondary = v;
        self
    }

    /// Exclude supplementary alignments (SAM flag 0x800).
    pub fn with_exclude_supplementary(mut self, v: bool) -> Self {
        self.exclude_supplementary = v;
        self
    }

    /// Set the strand to keep.
    pub fn with_strand(mut self, v: bio_types::strand::Strand) -> Self {
        self.strand = v;
        self
    }

    /// Keep only properly paired reads.
    pub fn with_proper_pair(mut self, v: bool) -> Self {
        self.proper_pair = v;
        self
    }

    /// Set the minimum mapping quality.
    pub fn with_min_mapq(mut self, v: u8) -> Self {
        self.min_mapq = v;
        self
    }

    /// Set the minimum read length.
    pub fn with_min_length(mut self, v: u32) -> Self {
        self.min_length = v;
        self
    }

    /// Set the maximum read length.
    pub fn with_max_length(mut self, v: u32) -> Self {
        self.max_length = v;
        self
    }

    /// Set the read group to keep.
    pub fn with_read_group(mut self, v: Option<String>) -> Self {
        self.read_group = v;
        self
    }

    /// Set the blacklisted locations (chromosome ID -> `Lapper`).
    pub fn with_blacklisted_locations(
        mut self,
        v: Option<HashMap<usize, Lapper<usize, u32>>>,
    ) -> Self {
        self.blacklisted_locations = v;
        self
    }

    /// Set the whitelisted cell barcodes.
    pub fn with_whitelisted_barcodes(mut self, v: Option<HashSet<String>>) -> Self {
        self.whitelisted_barcodes = v;
        self
    }

    /// Set the SAM tag to filter by.
    pub fn with_filter_tag(mut self, v: Option<String>) -> Self {
        self.filter_tag = v;
        self
    }

    /// Set the required value for `filter_tag`.
    pub fn with_filter_tag_value(mut self, v: Option<String>) -> Self {
        self.filter_tag_value = v;
        self
    }

    /// Set the minimum fragment length (template/insert size).
    pub fn with_min_fragment_length(mut self, v: Option<u32>) -> Self {
        self.min_fragment_length = v;
        self
    }

    /// Set the maximum fragment length (template/insert size).
    pub fn with_max_fragment_length(mut self, v: Option<u32>) -> Self {
        self.max_fragment_length = v;
        self
    }

    /// Checks if a read is valid according to the filter criteria.
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

        // Filter by duplicate/secondary/supplementary flags (cheapest checks first).
        if self.exclude_duplicates && flags.is_duplicate() {
            self.stats
                .n_failed_duplicate
                .fetch_add(1, Ordering::Relaxed);
            return Ok(false);
        }
        if self.exclude_secondary && flags.is_secondary() {
            self.stats
                .n_failed_secondary
                .fetch_add(1, Ordering::Relaxed);
            return Ok(false);
        }
        if self.exclude_supplementary && flags.is_supplementary() {
            self.stats
                .n_failed_supplementary
                .fetch_add(1, Ordering::Relaxed);
            return Ok(false);
        }

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
            self.stats
                .n_failed_proper_pair
                .fetch_add(1, Ordering::Relaxed);
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

        // Filter by fragment length (template length / insert size).
        if self.min_fragment_length.is_some() || self.max_fragment_length.is_some() {
            let tlen = alignment
                .template_length()
                .map_err(|e| anyhow::anyhow!(e))?
                .unsigned_abs();
            if let Some(min_flen) = self.min_fragment_length
                && tlen < min_flen
            {
                self.stats
                    .n_failed_fragment_length
                    .fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }
            if let Some(max_flen) = self.max_fragment_length
                && tlen > max_flen
            {
                self.stats
                    .n_failed_fragment_length
                    .fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }
        }

        let header = header.ok_or_else(|| anyhow::anyhow!("No header provided"))?;
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
        if let Some(blacklisted_locations) = &self.blacklisted_locations
            && let Some(blacklist) = blacklisted_locations.get(&chrom_id)
            && blacklist.count(start, end) > 0
        {
            self.stats
                .n_failed_blacklist
                .fetch_add(1, Ordering::Relaxed);
            return Ok(false);
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
            let read_group_value = match data.get(&rg) {
                None => {
                    self.stats
                        .n_not_in_read_group
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
                Some(v) => v?,
            };

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

        // Filter by tag.
        if let (Some(filter_tag), Some(filter_tag_value)) =
            (&self.filter_tag, &self.filter_tag_value)
        {
            if filter_tag.len() != 2 {
                self.stats
                    .n_failed_tag_filter
                    .fetch_add(1, Ordering::Relaxed);
                return Ok(false);
            }

            let tag_bytes = filter_tag.as_bytes();
            let tag = Tag::new(tag_bytes[0], tag_bytes[1]);
            let data = alignment.data();
            let tag_value = data.get(&tag);

            match tag_value {
                Some(Ok(Value::String(value))) => {
                    if value != filter_tag_value {
                        self.stats
                            .n_failed_tag_filter
                            .fetch_add(1, Ordering::Relaxed);
                        return Ok(false);
                    }
                }
                _ => {
                    self.stats
                        .n_failed_tag_filter
                        .fetch_add(1, Ordering::Relaxed);
                    return Ok(false);
                }
            }
        }

        Ok(true)
    }

    /// Returns a snapshot of the current filter statistics.
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

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::Flags;

    fn record_with_flags(flags: Flags) -> RecordBuf {
        let mut record = RecordBuf::default();
        *record.flags_mut() = flags;
        record
    }

    #[test]
    fn test_exclude_duplicates_filters_duplicate_reads() {
        let filter = BamReadFilter::default().with_exclude_duplicates(true);
        let record = record_with_flags(Flags::DUPLICATE);
        let valid = filter.is_valid(&record, None).expect("is_valid errored");
        assert!(!valid);
        assert_eq!(filter.stats().n_failed_duplicate(), 1);
    }

    #[test]
    fn test_exclude_secondary_filters_secondary_reads() {
        let filter = BamReadFilter::default().with_exclude_secondary(true);
        let record = record_with_flags(Flags::SECONDARY);
        let valid = filter.is_valid(&record, None).expect("is_valid errored");
        assert!(!valid);
        assert_eq!(filter.stats().n_failed_secondary(), 1);
    }

    #[test]
    fn test_exclude_supplementary_filters_supplementary_reads() {
        let filter = BamReadFilter::default().with_exclude_supplementary(true);
        let record = record_with_flags(Flags::SUPPLEMENTARY);
        let valid = filter.is_valid(&record, None).expect("is_valid errored");
        assert!(!valid);
        assert_eq!(filter.stats().n_failed_supplementary(), 1);
    }

    #[test]
    fn test_builder_methods_chain_and_apply() {
        let filter = BamReadFilter::default()
            .with_min_mapq(30)
            .with_min_length(50)
            .with_max_length(500)
            .with_proper_pair(false)
            .with_strand(bio_types::strand::Strand::Forward)
            .with_read_group(Some("RG1".to_string()))
            .with_filter_tag(Some("VP".to_string()))
            .with_filter_tag_value(Some("BCL2".to_string()))
            .with_min_fragment_length(Some(100))
            .with_max_fragment_length(Some(200))
            .with_exclude_duplicates(true);

        // Duplicate exclusion is checked first, so a duplicate read is rejected
        // regardless of the other builder settings above.
        let record = record_with_flags(Flags::DUPLICATE);
        let valid = filter.is_valid(&record, None).expect("is_valid errored");
        assert!(!valid);
        assert_eq!(filter.stats().n_failed_duplicate(), 1);
    }

    #[test]
    fn test_exclude_flags_default_off_does_not_filter() {
        let filter = BamReadFilter::default();
        let record = record_with_flags(Flags::DUPLICATE | Flags::SECONDARY | Flags::SUPPLEMENTARY);
        // The record is missing a header/alignment position, so the overall filter may
        // still fail or error further down the chain; we only care that the new
        // duplicate/secondary/supplementary counters stay untouched when disabled.
        let _ = filter.is_valid(&record, None);
        assert_eq!(filter.stats().n_failed_duplicate(), 0);
        assert_eq!(filter.stats().n_failed_secondary(), 0);
        assert_eq!(filter.stats().n_failed_supplementary(), 0);
    }
}
