use ahash::{HashMap, HashSet};
use noodles::{bam, sam::alignment::record::data::field::Value};
use rust_lapper::Lapper;
use std::fmt::Display;
use std::sync::{Arc, Mutex};
use anyhow::Result;

use crate::utils::CB;

#[derive(Clone, Copy, Debug)]
pub struct BamReadFilterStats {
    // Total number of reads
    n_total: u64,
    // Number of reads filtered by proper pair
    n_failed_proper_pair: u64,
    // Number of reads filtered by mapping quality
    n_failed_mapq: u64,
    // Number of reads filtered by length
    n_failed_length: u64,
    // Number of reads filtered by blacklisted locations
    n_failed_blacklist: u64,
    // Number of reads filtered by barcode
    n_failed_barcode: u64,
}

impl BamReadFilterStats {
    pub fn new() -> Self {
        Self {
            n_total: 0,
            n_failed_proper_pair: 0,
            n_failed_mapq: 0,
            n_failed_length: 0,
            n_failed_blacklist: 0,
            n_failed_barcode: 0,
        }
    }
}

/// A filter for BAM reads.
/// Set the minimum mapping quality, minimum and maximum read length, blacklisted locations, and whitelisted barcodes.
/// The filter is applied to each read in the BAM file.
#[derive(Debug)]
pub struct BamReadFilter {
    // Filter parameters

    // Properly paired reads only
    proper_pair: bool,
    // Minimum mapping quality
    min_mapq: u8,
    // Minimum read length
    min_length: u32,
    // Maximum read length
    max_length: u32,
    // Blacklisted locations (chromosome -> Lapper)
    blacklisted_locations: Option<HashMap<usize, Lapper<usize, i32>>>,
    // Whitelisted barcodes (cell barcodes to keep)
    whitelisted_barcodes: Option<HashSet<String>>,
    // Statistics for the filtering process
    // Arc<Mutex<>> is used to allow multiple threads to access the stats
    stats: Arc<Mutex<BamReadFilterStats>>,
}

impl Display for BamReadFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "\tProper pair: {}", self.proper_pair)?;
        writeln!(f, "\tMinimum mapping quality: {}", self.min_mapq)?;
        writeln!(f, "\tMinimum read length: {}", self.min_length)?;
        writeln!(f, "\tMaximum read length: {}", self.max_length)?;

        match &self.blacklisted_locations {
            Some(blacklisted_locations) => {
                writeln!(f, "\tNumber of Blacklisted locations: {:?}", blacklisted_locations.len())?;
            }
            None => {
                writeln!(f, "\tNumber of Blacklisted locations: 0")?;
            }
        }

        match &self.whitelisted_barcodes {
            Some(whitelisted_barcodes) => {
                writeln!(f, "\tNumber of Whitelisted barcodes: {:?}", whitelisted_barcodes.len())?;
            }
            None => {
                writeln!(f, "\tNumber of Whitelisted barcodes: 0")?;
            }
        }
        Ok(())
    }
}

impl Display for BamReadFilterStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "")?;
        writeln!(f, "Total reads: {}", self.n_total)?;
        writeln!(f, "Failed proper pair: {}", self.n_failed_proper_pair)?;
        writeln!(f, "Failed mapping quality: {}", self.n_failed_mapq)?;
        writeln!(f, "Failed length: {}", self.n_failed_length)?;
        writeln!(f, "Failed blacklist: {}", self.n_failed_blacklist)?;
        writeln!(f, "Failed barcode: {}", self.n_failed_barcode)?;
        Ok(())
    }
}


impl BamReadFilter {
    pub fn new(
        proper_pair: bool,
        min_mapq: Option<u8>,
        min_length: Option<u32>,
        max_length: Option<u32>,
        blacklisted_locations: Option<HashMap<usize, Lapper<usize, i32>>>,
        whitelisted_barcodes: Option<HashSet<String>>,
    ) -> Self {
        let min_mapq = min_mapq.unwrap_or(0);
        let min_length = min_length.unwrap_or(0);
        let max_length = max_length.unwrap_or(std::u32::MAX);

        Self {
            proper_pair,
            min_mapq,
            min_length,
            max_length,
            blacklisted_locations,
            whitelisted_barcodes,
            stats: Arc::new(Mutex::new(BamReadFilterStats::new())),
        }
    }

    pub fn is_valid(&self, alignment: &bam::Record) -> bool {
        // Update total
        self.stats.lock().unwrap().n_total += 1;


        // Filter by proper pair
        if self.proper_pair && !alignment.flags().is_properly_segmented() {
            return false;
        }

        // Filter by unmapped reads
        if alignment.flags().is_unmapped() {
            self.stats.lock().unwrap().n_failed_mapq += 1;
            return false;
        }

        // Filter by mapping quality
        if alignment.mapping_quality().unwrap().get() < self.min_mapq {
            self.stats.lock().unwrap().n_failed_mapq += 1;
            return false;
        }

        // Filter by read length
        let alignment_length = alignment.sequence().len();

        if alignment_length < self.min_length as usize
            || alignment_length > self.max_length as usize
        {
            self.stats.lock().unwrap().n_failed_length += 1;
            return false;
        }

        let chrom = alignment
            .reference_sequence_id()
            .expect("No reference sequence ID")
            .expect("Failed to get reference sequence ID");
        let start = alignment
            .alignment_start()
            .expect("No alignment start")
            .expect("Failed to get alignment start")
            .get();
        let end = start + alignment_length;

        // Filter by blacklisted locations
        match &self.blacklisted_locations {
            Some(blacklisted_locations) => {
                if let Some(blacklist) = blacklisted_locations.get(&chrom) {
                    if blacklist.count(start, end) > 0 {
                        self.stats.lock().unwrap().n_failed_blacklist += 1;
                        return false;
                    }
                }
            }
            None => {}
        }

        // Filter by blacklisted barcodes
        if let Some(barcodes) = &self.whitelisted_barcodes {
            let barcode = get_cell_barcode(alignment);

            match barcode {
                Some(barcode) => {
                    if !barcodes.contains(&barcode) {
                        self.stats.lock().unwrap().n_failed_barcode += 1;
                        return false;
                    }
                }
                None => {
                    self.stats.lock().unwrap().n_failed_barcode += 1;
                    return false;
                }
            }
        };
        true
    }

    pub fn stats(&self) -> BamReadFilterStats {
        let stats = self.stats.lock().unwrap();
        *stats
    }


}

/// Get the cell barcode from a BAM alignment.
fn get_cell_barcode(alignment: &bam::Record) -> Option<String> {
    let tags = alignment.data();
    let cb = tags.get(&CB);

    match cb {
        Some(Ok(Value::String(barcode))) => Some(barcode.to_string()),
        _ => None,
    }

}
