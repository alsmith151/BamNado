use crate::read_filter::BamReadFilter;
use crate::bam_utils::get_bam_header;
use crate::bam_utils::BamStats;
use anyhow::Result;
use crossbeam::channel::unbounded;
use log::{debug, info, warn};
use indicatif::ProgressBar;


use noodles::core::{Position, Region};
use noodles::{bam, sam};
use noodles::bam::bai;
use noodles::sam::header::record::value::{map::ReferenceSequence, Map};
use noodles::bam::r#async::io::{Reader as AsyncReader, Writer as AsyncWriter};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::path::PathBuf;

use futures::TryStreamExt;
use tokio::fs::File;

pub struct BamFilterer {
    filepath: PathBuf,
    output: PathBuf,
    filter: BamReadFilter,
}

impl BamFilterer {
    pub fn new(filepath: PathBuf, output: PathBuf, filter: BamReadFilter) -> Self {
        BamFilterer {
            filepath,
            output,
            filter,
        }
    }

    pub async fn split_async(&self) -> Result<()> {
        let filepath = self.filepath.clone();
        let outfile = self.output.clone();

        let mut reader = File::open(&filepath).await.map(AsyncReader::new)?;
        let header = reader.read_header().await;

        let header = match header {
            Ok(header) => header,
            Err(e) => {
               get_bam_header(&filepath)?
            }
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
                if is_valid {
                    writer.write_record(&header, &record).await?;
                }
            }
        }

        progress.finish();
        writer.shutdown().await?;



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
                    .filter(|record| {
                        self.filter
                            .is_valid(record, Some(&header))
                            .unwrap_or_else(|_| false)
                    })
                    .collect::<Vec<_>>();

                tx.send(filtered_records).expect("Error sending records");
            });

        handle.join().expect("Error joining threads");
        Ok(())
    }
}
