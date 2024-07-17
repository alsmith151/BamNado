use anyhow::Ok;
use bio::io::bed;
use log::{info, warn};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Format, Header, IndexedReader, Read};
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
type Iv = Interval<u64, u64>;
use std::str;
use std::sync::Arc;
use std::thread;

fn get_intervals(bed: &PathBuf) -> Result<HashMap<String, Vec<Iv>>, anyhow::Error> {
    let mut bed_intervals = HashMap::new();
    let mut reader = bed::Reader::from_file(Path::new(&bed)).expect("Could not open BED file");

    for record in reader.records() {
        let record = record.expect("Error reading BED record");
        let interval = Iv {
            start: record.start(),
            stop: record.end(),
            val: 0,
        };

        let chrom = record.chrom().to_owned();
        bed_intervals.entry(chrom).or_insert(vec![]).push(interval);
    }

    Ok(bed_intervals)
}

fn get_chrom_names(header: &rust_htslib::bam::HeaderView) -> Result<Vec<String>, anyhow::Error> {
    let tids: Vec<_> = header
        .target_names()
        .iter()
        .map(|n| header.tid(n))
        .filter(|tid| tid.is_some())
        .map(|tid| {
            str::from_utf8(header.tid2name(tid.unwrap()))
                .expect("Could not convert chrom name to str")
                .to_owned()
        })
        .collect();
    Ok(tids)
}

pub fn remove_regions_from_bam(
    bed: PathBuf,
    bam: PathBuf,
    output: PathBuf,
    n_threads: usize,
) -> Result<(), anyhow::Error> {
    let intervals_for_subtraction =
        Arc::new(get_intervals(&bed).expect("Could not get intervals from BED file"));

    let bam_reader = rust_htslib::bam::Reader::from_path(&bam).expect("Could not open BAM file");
    let header_view = bam_reader.header().to_owned();
    let header = Header::from_template(&header_view);
    let chrom_names = get_chrom_names(&header_view).expect("Could not get chrom names");

    let (chrom_sender, chrom_recv) = crossbeam::channel::unbounded::<String>();
    let (filt_sender, filt_recv) = crossbeam::channel::unbounded();

    let mut filter_handles = Vec::new();

    // Spawn filtering threads
    for _ in 0..n_threads {
        let chrom_recv = chrom_recv.clone();
        let writer_sender = filt_sender.clone();
        let intervals_for_subtraction = intervals_for_subtraction.clone();
        let bam = bam.clone();

        filter_handles.push(thread::spawn(move || {
            for chrom in chrom_recv {
                let mut record_batch = Vec::with_capacity(1e5 as usize);
                let mut batch_counter = 0;

                match intervals_for_subtraction.get(&chrom) {
                    Some(intervals) => {
                        let mut reader =
                            IndexedReader::from_path(&bam).expect("Could not open BAM file");
                        reader.fetch(&chrom).expect("Failed to fetch chromosome");

                        let lapper = Lapper::new(intervals.clone());

                        for result in reader.records() {
                            if batch_counter == 1e5 as usize {
                                writer_sender
                                    .send(record_batch)
                                    .expect("Failed to send records");
                                record_batch = Vec::with_capacity(1e5 as usize);
                                batch_counter = 0;
                            }

                            let record = result.expect("Could not read BAM record");
                            let start = record.reference_start() as u64;
                            let end = record.reference_end() as u64;

                            let overlap_count = lapper.count(start, end);
                            if overlap_count == 0 {
                                record_batch.push(record);
                                batch_counter += 1;
                            }
                        }

                        // Send any remaining records
                        if record_batch.len() > 0 {
                            writer_sender
                                .send(record_batch)
                                .expect("Failed to send records");
                        }
                    }
                    None => {
                        let mut reader =
                            IndexedReader::from_path(&bam).expect("Could not open BAM file");
                        reader.fetch(&chrom).expect("Failed to fetch chromosome");

                        for result in reader.records() {
                            if batch_counter == 1e5 as usize {
                                writer_sender
                                    .send(record_batch)
                                    .expect("Failed to send records");
                                record_batch = Vec::with_capacity(1e5 as usize);
                                batch_counter = 0;
                            }

                            let record = result.expect("Could not read BAM record");
                            record_batch.push(record);
                            batch_counter += 1;
                        }

                        // Send any remaining records
                        if record_batch.len() > 0 {
                            writer_sender
                                .send(record_batch)
                                .expect("Failed to send records");
                        }
                    }
                }
            }
            // Drop the sender so the receiver will know we're done
            drop(writer_sender);
        }));
    }

    // Spawn writing thread
    let writer_handle = thread::spawn(move || {
        let mut bam_writer = rust_htslib::bam::Writer::from_path(output, &header, Format::Bam)
            .expect("Could not open BAM file for writing");

        for record_batch in filt_recv {
            for read in record_batch {
                bam_writer.write(&read).expect("Failed to write record");
            }
        }
    });

    // Send chromosomes to threads
    for chrom in chrom_names {
        chrom_sender.send(chrom)?;
    }

    // Drop the sender so the receiver will know we're done
    drop(chrom_sender);
    drop(filt_sender); // Drop the ref to the sender so the threads will know we're done

    // Join threads
    for handle in filter_handles {
        handle.join().expect("Failed to join filter thread");
    }
    writer_handle.join().expect("Failed to join writer thread");

    Ok(())
}

// Test remove regions from bam
#[cfg(test)]
#[test]
fn test_remove_regions_from_bam() {
    let bed = PathBuf::from("test/test_subtraction.bed");
    let bam = PathBuf::from("test/iALL-863388_H3K27ac-1_subsample.bam");
    let output = PathBuf::from("test/test_no_regions.bam");
    let n_threads = 4;

    remove_regions_from_bam(bed, bam, output, n_threads)
        .expect("Could not remove regions from BAM file");
}
