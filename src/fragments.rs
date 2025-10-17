use ahash::{HashMap, HashMapExt, HashSet};
use anyhow::{Context, Result};
use noodles::sam::alignment::record::data::field::tag::Tag;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::Record;
use regex::Regex;
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Configuration for ATAC-seq fragment extraction
#[derive(Debug, Clone)]
pub struct FragmentConfig {
    /// Minimum mapping quality for reads
    pub min_mapq: u8,
    /// Cell barcode tag (e.g., "CB")
    pub cell_barcode_tag: Tag,
    /// Optional regex to extract barcode from read name
    pub readname_barcode_regex: Option<Regex>,
    /// Optional list of cell barcodes to retain
    pub cells: Option<HashSet<String>>,
    /// Maximum distance between fragment ends (bp)
    pub max_distance: i32,
    /// Minimum distance between fragment ends (bp)
    pub min_distance: i32,
    /// Number of reads to process before collapsing fragments
    pub chunk_size: usize,
    /// Position shifts for Tn5 correction: [forward_strand, reverse_strand]
    pub shifts: [i32; 2],
    /// Only collapse fragments within the same cell barcode
    pub collapse_within_barcode: bool,
    /// Maximum distance for collapsing overlapping fragments
    pub max_collapse_distance: i32,
}

impl Default for FragmentConfig {
    fn default() -> Self {
        Self {
            min_mapq: 30,
            cell_barcode_tag: Tag::from([b'C', b'B']), // CB tag
            readname_barcode_regex: None,
            cells: None,
            max_distance: 5000,
            min_distance: 10,
            chunk_size: 500_000,
            shifts: [4, -5],
            collapse_within_barcode: false,
            max_collapse_distance: 20,
        }
    }
}

/// Represents a single ATAC-seq fragment
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Fragment {
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub barcode: String,
    pub count: u32,
}

impl Fragment {
    /// Create a unique key for fragment without barcode (for deduplication across cells)
    fn coord_key(&self) -> String {
        format!("{}|{}|{}", self.chromosome, self.start, self.end)
    }

    /// Create a unique key for fragment with barcode (for deduplication within cells)
    fn full_key(&self) -> String {
        format!(
            "{}|{}|{}|{}",
            self.chromosome, self.start, self.end, self.barcode
        )
    }

    /// Get key with start position only (for overlap detection)
    fn start_key(&self) -> String {
        format!("{}|{}|{}", self.chromosome, self.start, self.barcode)
    }

    /// Get key with end position only (for overlap detection)
    fn end_key(&self) -> String {
        format!("{}|{}|{}", self.chromosome, self.end, self.barcode)
    }

    /// Write fragment to output
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            self.chromosome, self.start, self.end, self.barcode, self.count
        )
        .context("Failed to write fragment")
    }
}

/// Represents an incomplete fragment being assembled from paired-end reads
#[derive(Debug, Clone)]
struct IncompleteFragment {
    chromosome: String,
    start: Option<i32>,
    end: Option<i32>,
    barcode: Option<String>,
}

impl IncompleteFragment {
    fn new(chromosome: String) -> Self {
        Self {
            chromosome,
            start: None,
            end: None,
            barcode: None,
        }
    }

    /// Check if fragment is complete (has both start and end)
    fn is_complete(&self) -> bool {
        self.start.is_some() && self.end.is_some()
    }

    /// Get the furthest position in the fragment (for determining completion)
    fn max_position(&self) -> Option<i32> {
        match (self.start, self.end) {
            (Some(s), Some(e)) => Some(max(s, e)),
            (Some(s), None) => Some(s),
            (None, Some(e)) => Some(e),
            (None, None) => None,
        }
    }

    /// Convert to a complete Fragment
    fn into_fragment(self) -> Option<Fragment> {
        match (self.start, self.end, self.barcode) {
            (Some(start), Some(end), Some(barcode)) => Some(Fragment {
                chromosome: self.chromosome,
                start: min(start, end),
                end: max(start, end),
                barcode,
                count: 1,
            }),
            _ => None,
        }
    }
}

/// Manager for tracking incomplete fragments during BAM processing
pub struct FragmentTracker {
    /// Map of read name to incomplete fragment
    incomplete: HashMap<String, IncompleteFragment>,
    /// Completed fragments ready for deduplication
    complete: Vec<Fragment>,
    /// Configuration
    config: FragmentConfig,
    /// Statistics
    reads_processed: usize,
    fragments_created: usize,
    fragments_filtered: usize,
}

impl FragmentTracker {
    pub fn new(config: FragmentConfig) -> Self {
        Self {
            incomplete: HashMap::new(),
            complete: Vec::new(),
            config,
            reads_processed: 0,
            fragments_created: 0,
            fragments_filtered: 0,
        }
    }

    /// Process a single alignment record
    pub fn process_record<R: Record>(&mut self, record: &R, chromosome: &String) -> Result<()> {
        self.reads_processed += 1;

        // Extract mapping quality
        let mapq = match record.mapping_quality() {
            Some(Ok(mq)) => u8::from(mq),
            _ => 0,
        };
        
        if mapq < self.config.min_mapq {
            return Ok(());
        }

        // Extract positions
        let alignment_start = match record.alignment_start() {
            Some(Ok(pos)) => usize::from(pos) as i32 - 1, // Convert to 0-based
            _ => return Ok(()),
        };
        
        let alignment_end = match record.alignment_end() {
            Some(Ok(pos)) => usize::from(pos) as i32,
            _ => return Ok(()),
        };

        // Get flags and check strand
        let is_reverse = record.flags()?.is_reverse_complemented();
        let (start, end) = if is_reverse {
            (alignment_start, alignment_end + self.config.shifts[1])
        } else {
            (alignment_start + self.config.shifts[0], alignment_end)
        };

        // Extract cell barcode
        let barcode = self.extract_barcode(record)?;

        // Filter by cell barcode if specified
        if let Some(ref allowed_cells) = self.config.cells {
            if let Some(ref bc) = barcode {
                if !allowed_cells.contains(bc) {
                    return Ok(());
                }
            } else {
                return Ok(());
            }
        }

        // Update fragment tracker
        let read_name = std::str::from_utf8(record.name().unwrap_or_default())
            .unwrap_or("")
            .to_string();

        self.update_fragment(read_name, chromosome.to_string(), start, end, record.flags()?.is_first_segment(), record.flags()?.is_reverse_complemented(), barcode)?;

        Ok(())
    }


    /// Extract cell barcode from record
    fn extract_barcode<R: Record>(&self, record: &R) -> Result<Option<String>> {
        // Try regex extraction from read name first
        if let Some(ref regex) = self.config.readname_barcode_regex {
            if let Some(name_bytes) = record.name() {
                let name = std::str::from_utf8(name_bytes)?;
                if let Some(captures) = regex.captures(name) {
                    if let Some(matched) = captures.get(0) {
                        return Ok(Some(matched.as_str().to_string()));
                    }
                }
            }
        }

        // Otherwise, use cell barcode tag
        let data = record.data();
        if let Some(Ok(Value::String(bc))) = data.get(&self.config.cell_barcode_tag) {
            return Ok(Some(bc.to_string()));
        }

        Ok(None)
    }

    /// Update fragment information with a new read
    fn update_fragment(
        &mut self,
        read_name: String,
        chromosome: String,
        start: i32,
        end: i32,
        is_first_segment: bool,
        is_reverse_complemented: bool,
        barcode: Option<String>,
    ) -> Result<()> {
        
        

        if let Some(mut frag) = self.incomplete.remove(&read_name) {
            // Fragment already exists, add the mate information
            match (is_first_segment, is_reverse_complemented) {
                // First read of the pair and forward strand
                (true, false) => {
                    frag.start = Some(start);
                }
                // First read of the pair and reverse strand
                (true, true) => {
                    frag.end = Some(end);
                }
                // Second read of the pair and forward strand
                (false, false) => {
                    frag.end = Some(end);
                }
                // Second read of the pair and reverse strand
                (false, true) => {
                    frag.start = Some(start);
                }
            }
                
            // Update barcode if we have one and don't already
            if frag.barcode.is_none() && barcode.is_some() {
                frag.barcode = barcode;
            }

            // Check if fragment is complete and valid
            if frag.is_complete() {
                if let Some(fragment) = frag.into_fragment() {
                    let fragment_len = fragment.end - fragment.start;
                    
                    // Validate fragment size
                    if fragment_len >= self.config.min_distance
                        && fragment_len <= self.config.max_distance
                    {
                        self.complete.push(fragment);
                        self.fragments_created += 1;
                    } else {
                        self.fragments_filtered += 1;
                    }
                } else {
                    // Missing barcode or other issue
                    self.fragments_filtered += 1;
                }
            } else {
                // Put back incomplete fragment
                self.incomplete.insert(read_name, frag);
            }
        } else {

            // New fragment - this is the first read of the pair
            let mut frag = IncompleteFragment::new(chromosome);
            frag.start = Some(start);
            frag.end = Some(end); 
            frag.barcode = barcode;
            self.incomplete.insert(read_name, frag);
        }

        Ok(())
    }

    /// Remove and return completed fragments that are far enough from current position
    pub fn harvest_complete_fragments(&mut self, current_position: i32) -> Vec<Fragment> {
        let threshold = current_position - self.config.max_distance - self.config.max_collapse_distance;

        // Split complete fragments into those ready to harvest and those to keep
        let (ready, keep): (Vec<Fragment>, Vec<Fragment>) = self.complete
            .drain(..)
            .partition(|frag| frag.end < threshold);
        
        self.complete = keep;

        // Remove incomplete fragments that are too old
        self.incomplete.retain(|_, frag| {
            if let Some(max_pos) = frag.max_position() {
                max_pos >= threshold
            } else {
                false // Remove fragments with no position set
            }
        });

        ready
    }

    /// Get all remaining fragments (call at end of chromosome)
    pub fn harvest_all(&mut self) -> Vec<Fragment> {
        // Convert any remaining incomplete fragments with both start and end
        let remaining: Vec<_> = self
            .incomplete
            .drain()
            .filter_map(|(_, frag)| frag.into_fragment())
            .collect();

        let mut all_fragments = std::mem::take(&mut self.complete);
        all_fragments.extend(remaining);
        all_fragments
    }

    pub fn stats(&self) -> (usize, usize, usize) {
        (
            self.reads_processed,
            self.fragments_created,
            self.fragments_filtered,
        )
    }
}

/// Collapse duplicate and overlapping fragments
pub struct FragmentCollapser {
    config: FragmentConfig,
}

impl FragmentCollapser {
    pub fn new(config: FragmentConfig) -> Self {
        Self { config }
    }

    /// Collapse fragments with deduplication and overlap handling
    pub fn collapse(&self, fragments: Vec<Fragment>) -> Vec<Fragment> {
        if fragments.is_empty() {
            return Vec::new();
        }

        // Group fragments by full key and count duplicates
        let mut fragment_counts: HashMap<String, Fragment> = HashMap::new();
        for frag in fragments {
            let key = if self.config.collapse_within_barcode {
                frag.full_key()
            } else {
                frag.coord_key()
            };

            fragment_counts
                .entry(key)
                .and_modify(|f| f.count += 1)
                .or_insert(frag);
        }

        let mut fragments: Vec<_> = fragment_counts.into_values().collect();

        // If not collapsing within barcodes, handle overlapping fragments
        if !self.config.collapse_within_barcode {
            fragments = self.collapse_overlaps(fragments);
        }

        fragments
    }

    /// Collapse fragments that share start or end coordinates within the same barcode
    fn collapse_overlaps(&self, mut fragments: Vec<Fragment>) -> Vec<Fragment> {
        // Sort by chromosome, then position
        fragments.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });

        // Collapse by shared start positions
        fragments = self.collapse_by_position(fragments, true);
        
        // Collapse by shared end positions
        fragments = self.collapse_by_position(fragments, false);

        fragments
    }

    /// Collapse fragments sharing a start (pos=true) or end (pos=false) coordinate
    fn collapse_by_position(&self, fragments: Vec<Fragment>, use_start: bool) -> Vec<Fragment> {
        // Group by position key
        let mut position_groups: HashMap<String, Vec<Fragment>> = HashMap::new();
        
        for frag in fragments {
            let key = if use_start {
                frag.start_key()
            } else {
                frag.end_key()
            };
            position_groups.entry(key).or_default().push(frag);
        }

        let mut result = Vec::new();

        for (_, mut group) in position_groups {
            if group.len() == 1 {
                result.push(group.into_iter().next().unwrap());
            } else {
                // Find fragment with maximum count
                group.sort_by_key(|f| std::cmp::Reverse(f.count));
                let mut winner = group[0].clone();
                
                // Sum all counts
                let total_count: u32 = group.iter().map(|f| f.count).sum();
                winner.count = total_count;
                
                result.push(winner);
            }
        }

        result
    }
}



/// Write fragments to file
pub fn write_fragments<P: AsRef<Path>>(
    fragments: &[Fragment],
    path: P,
    append: bool,
) -> Result<()> {
    
    
    let file = if append {
        File::options()
            .create(true)
            .append(true)
            .open(&path)?
    } else {
        File::create(&path)?
    };

    // Change the writer to support compression
    let mut writer: Box<dyn Write> = match path.as_ref().extension().and_then(|s| s.to_str()) {
        Some("zst") => {
            let encoder = zstd::stream::Encoder::new(file, 0)?;
            Box::new(BufWriter::new(encoder))
        },
        Some("gz") => {
            let encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            Box::new(BufWriter::new(encoder))
        },
        _ => Box::new(BufWriter::new(file)),
    };
   
    for fragment in fragments {
        fragment.write(&mut writer)?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fragment_keys() {
        let frag = Fragment {
            chromosome: "chr1".to_string(),
            start: 1000,
            end: 1200,
            barcode: "ACGT".to_string(),
            count: 1,
        };

        assert_eq!(frag.coord_key(), "chr1|1000|1200");
        assert_eq!(frag.full_key(), "chr1|1000|1200|ACGT");
        assert_eq!(frag.start_key(), "chr1|1000|ACGT");
        assert_eq!(frag.end_key(), "chr1|1200|ACGT");
    }

    #[test]
    fn test_fragment_collapse() {
        let config = FragmentConfig::default();
        let collapser = FragmentCollapser::new(config);

        let fragments = vec![
            Fragment {
                chromosome: "chr1".to_string(),
                start: 1000,
                end: 1200,
                barcode: "ACGT".to_string(),
                count: 1,
            },
            Fragment {
                chromosome: "chr1".to_string(),
                start: 1000,
                end: 1200,
                barcode: "ACGT".to_string(),
                count: 1,
            },
            Fragment {
                chromosome: "chr1".to_string(),
                start: 1000,
                end: 1200,
                barcode: "TGCA".to_string(),
                count: 1,
            },
        ];

        let collapsed = collapser.collapse(fragments);
        
        // Should collapse to one fragment with count=3 (same coordinates)
        assert_eq!(collapsed.len(), 1);
        assert_eq!(collapsed[0].count, 3);
    }

    #[test]
    fn test_incomplete_fragment() {
        let mut frag = IncompleteFragment::new("chr1".to_string());
        assert!(!frag.is_complete());

        frag.start = Some(1000);
        assert!(!frag.is_complete());

        frag.end = Some(1200);
        frag.barcode = Some("ACGT".to_string());
        assert!(frag.is_complete());

        let complete = frag.into_fragment().unwrap();
        assert_eq!(complete.start, 1000);
        assert_eq!(complete.end, 1200);
    }
}

