use crate::filter::BamReadFilter;
use crate::utils::get_bam_header;
use crate::utils::BamStats;
use anyhow::Result;
use crossbeam::channel::unbounded;
use noodles::bam;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::path::PathBuf;

pub struct BamSplitter {
    filepath: PathBuf,
    output: PathBuf,
    filter: BamReadFilter,
}

impl BamSplitter {
    pub fn new(filepath: PathBuf, output: PathBuf, filter: BamReadFilter) -> Self {
        BamSplitter {
            filepath,
            output,
            filter,
        }
    }

    pub fn split(&self) -> Result<()> {
        let filepath = self.filepath.clone();

        // Determine the chunks over which to iterate
        // Assuming here that the depth of the first file is the same as the rest
        let bam_stats = BamStats::new(self.filepath.clone())?;
        let genomic_chunks = bam_stats.genome_chunks(1e6 as u64)?;

        // Set up required variables
        let chromsizes_refid = bam_stats
            .chromsizes_ref_id()
            .expect("Error getting chromsizes");

        // Iterate over the genomic chunks and pileup the reads
        let n_total_chunks = genomic_chunks.len();

        // Start a writing thread
        let outfile = self.output.clone();

        // Set-up for multithreaded writing
        let (tx, rx) = unbounded();

        // Start a writing thread
        let handle = std::thread::spawn(move || {
            let reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(&filepath)
                .expect("Error reading file");
            let header = get_bam_header(&filepath).expect("Error reading header");
            let mut writer = bam::io::writer::Builder::default()
                .build_from_path(outfile)
                .expect("Error writing to file");

            for chunk in rx.iter() {
                for record in chunk {
                    writer
                        .write_record(&header, &record)
                        .expect("Error writing record");
                }
            }
        });

        let filepath = self.filepath.clone();

        // Process the genomic chunks filtering reads
        genomic_chunks
            .into_par_iter()
            .for_each_with(tx, |tx, chunk| {
                
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(&filepath)
                    .expect("Error reading file");
                
                let header = get_bam_header(&filepath).expect("Error reading header");
                
                let records = reader
                    .query(&header, &chunk)
                    .expect("Error getting chunk of reads");

                let filtered_records = records
                    .into_iter()
                    .filter_map(|r| r.is_ok().then(|| r.unwrap()))
                    .filter(|record| self.filter.is_valid(record))
                    .collect::<Vec<_>>();

                tx.send(filtered_records).expect("Error sending records");
            });

        handle.join().expect("Error joining threads");
        Ok(())
    }
}
