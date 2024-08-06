use crate::filter::BamReadFilter;
use crate::utils::get_bam_header;
use crate::utils::BamStats;
use anyhow::Result;
use crossbeam::channel::unbounded;
use itertools::Itertools;
use noodles::core::{Region, Position};
use noodles::{bam};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::path::PathBuf;
use log::{info, warn, debug};

use futures::TryStreamExt;
use tokio::fs::File;

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

    pub async fn split_async(&self) -> Result<()> {
        let filepath = self.filepath.clone();
        let outfile = self.output.clone();

        let mut reader = File::open(&filepath).await.map(bam::AsyncReader::new)?;
        let header = reader.read_header().await;

        let header = match header {
            Ok(header) => header,
            Err(e) => get_bam_header(filepath.clone())?,
        }; 

        let index = bam::bai::r#async::read(self.filepath.with_extension("bam.bai")).await?;

        // Make writer
        let mut writer = bam::AsyncWriter::new(File::create(outfile).await?);
        writer.write_header(&header).await?;

        // Get the chromosome sizes
        let chromsizes = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)|
             (name.to_string(), seq.length().get() as u64)
            )
            .collect::<std::collections::HashMap<_, _>>();

        let query_regions = chromsizes
            .iter()
            .map(|(name, size)| {
                let start = Position::try_from(1).unwrap();
                let end = Position::try_from(*size as usize).unwrap();
                Region::new(name.to_string(), start..=end)

            });


        for region in query_regions {
            let mut query = reader.query(&header, &index, &region)?;
            while let Some(record) =  query.try_next().await? {
                let is_valid = self.filter.is_valid(&record, Some(&header))?;
                if is_valid {
                    writer.write_record(&header, &record).await?;
                }

                
            }
        }




        // loop  {

        //     let record = records.try_next().await;

        //     // Check if the record is valid i.e. doesn't error
        //     let rec = match record {
        //         Ok(rec) => rec,
        //         Err(e) => {
        //             warn!("Error reading record: {:?}", e);
        //             continue;
        //         }
        //     };

        //     // Check if the record is None
        //     let record = match rec {
        //         Some(record) => record,
        //         None => break,
                
        //     };

        //     // Check if the record is valid i.e. passes the filter and write it to the output file
        //     let is_valid = self.filter.is_valid(&record, Some(&header))?;
        //     if is_valid {
        //         writer.write_alignment_record(&header, &record).await?;
        //     }
        // }

        Ok(())
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
                    .filter(|record| self.filter.is_valid(record, Some(&header)).unwrap_or_else(|_| false))
                    .collect::<Vec<_>>();

                tx.send(filtered_records).expect("Error sending records");
            });

        handle.join().expect("Error joining threads");
        Ok(())
    }
}
