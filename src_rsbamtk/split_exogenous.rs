use ahash::HashMap;
use anyhow::{Context, Result};
use bstr::ByteSlice;
use clap::builder::Str;
use crossbeam::channel::{bounded, unbounded, Receiver, Sender};
use crossbeam::queue::ArrayQueue;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator};
use noodles::bam::bai::read;
use noodles::bam::io::Writer;
use noodles::bed::record;
use noodles::core::region;
use noodles::sam::header;
use noodles::{bam, bgzf, sam};
use rayon::prelude::*;
use sam::header::record::value::{map::ReferenceSequence, Map};
use serde::{Deserialize, Serialize};
use std::fmt::format;
use std::num::NonZeroUsize;
use std::ops::Add;
use std::path::PathBuf;
use std::prelude::v1::*;

#[derive(Debug, Serialize, Deserialize)]
pub struct SplitStats {
    filename: String,
    n_unmapped_reads: u64,
    n_qcfail_reads: u64,
    n_duplicate_reads: u64,
    n_secondary_reads: u64,
    n_low_maq: u64,
    n_both_genomes: u64,
    n_exogenous: u64,
    n_endogenous: u64,
}

impl SplitStats {
    fn new(filename: String) -> Self {
        Self {
            filename,
            n_unmapped_reads: 0,
            n_qcfail_reads: 0,
            n_duplicate_reads: 0,
            n_secondary_reads: 0,
            n_low_maq: 0,
            n_both_genomes: 0,
            n_exogenous: 0,
            n_endogenous: 0,
        }
    }

    fn add_unmapped(&mut self) {
        self.n_unmapped_reads += 1;
    }

    fn add_qcfail(&mut self) {
        self.n_qcfail_reads += 1;
    }

    fn add_duplicate(&mut self) {
        self.n_duplicate_reads += 1;
    }

    fn add_secondary(&mut self) {
        self.n_secondary_reads += 1;
    }

    fn add_low_maq(&mut self) {
        self.n_low_maq += 1;
    }

    fn add_both_genomes(&mut self) {
        self.n_both_genomes += 1;
    }

    fn add_exogenous(&mut self) {
        self.n_exogenous += 1;
    }

    fn add_endogenous(&mut self) {
        self.n_endogenous += 1;
    }

    pub fn print(&self) {
        println!("Filename: {}", self.filename);
        println!("Unmapped reads: {}", self.n_unmapped_reads);
        println!("QC fail reads: {}", self.n_qcfail_reads);
        println!("Duplicate reads: {}", self.n_duplicate_reads);
        println!("Secondary reads: {}", self.n_secondary_reads);
        println!("Low mapping quality reads: {}", self.n_low_maq);
        println!("Both genomes reads: {}", self.n_both_genomes);
        println!("Exogenous reads: {}", self.n_exogenous);
        println!("Endogenous reads: {}", self.n_endogenous);
    }
}

impl Add for SplitStats {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            filename: self.filename,
            n_unmapped_reads: self.n_unmapped_reads + other.n_unmapped_reads,
            n_qcfail_reads: self.n_qcfail_reads + other.n_qcfail_reads,
            n_duplicate_reads: self.n_duplicate_reads + other.n_duplicate_reads,
            n_secondary_reads: self.n_secondary_reads + other.n_secondary_reads,
            n_low_maq: self.n_low_maq + other.n_low_maq,
            n_both_genomes: self.n_both_genomes + other.n_both_genomes,
            n_exogenous: self.n_exogenous + other.n_exogenous,
            n_endogenous: self.n_endogenous + other.n_endogenous,
        }
    }
}

pub struct SplitBam {
    bam_input: PathBuf,
    bam_endogenous: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_exogenous: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_both_genomes: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
    bam_unmapped: bam::io::Writer<noodles::bgzf::Writer<std::fs::File>>,
}

struct BamHeaders {
    header_input: sam::Header,
    header_endogenous: sam::Header,
    header_exogenous: sam::Header,
    header_both_genomes: sam::Header,
    header_unmapped: sam::Header,
}

impl SplitBam {
    pub fn new(bam_input: PathBuf, output_prefix: PathBuf) -> Result<Self> {
        let bam_input = bam_input;
        let bam_endogenous = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("endogenous.bam"))?;
        let bam_exogenous = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("exogenous.bam"))?;
        let bam_both_genomes = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("both_genomes.bam"))?;
        let bam_unmapped = bam::io::writer::Builder::default()
            .build_from_path(output_prefix.with_extension("unmapped.bam"))?;

        Ok(Self {
            bam_input,
            bam_endogenous,
            bam_exogenous,
            bam_both_genomes,
            bam_unmapped,
        })
    }

    fn get_input_reader(&self) -> Result<bam::io::IndexedReader<bgzf::Reader<std::fs::File>>> {
        // Check if the input BAM file has an index file.
        let index_path = self.bam_input.with_extension("bam.bai");
        if !index_path.exists() {
            return Err(anyhow::anyhow!("Index file not found for BAM file"));
        }
        let reader = bam::io::indexed_reader::Builder::default().build_from_path(self.bam_input.clone())?;
        Ok(reader)
    }

    fn make_headers(&mut self, exogenous_prefix: &[u8]) -> Result<BamHeaders> {

        let header_input = self.get_input_reader().expect("Error reading bam input file").read_header()?;
        let reference_seqs = header_input.reference_sequences().clone();

        // Split reference sequences into endogenous and exogenous based on prefixes present.
        // Endogenous sequences have no prefix, exogenous sequences have a prefix.
        let mut reference_seqs_endogenous = sam::header::ReferenceSequences::new();
        let mut reference_seqs_exogenous = sam::header::ReferenceSequences::new();

        for (name, len) in reference_seqs.iter() {
            if name.starts_with(&exogenous_prefix) {
                reference_seqs_exogenous.insert(name.clone(), len.clone());
            } else {
                reference_seqs_endogenous.insert(name.clone(), len.clone());
            }
        }

        let header_endogenous = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs_endogenous)
            .build();

        let header_exogenous = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs_exogenous)
            .build();

        let header_both_genomes = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs.clone())
            .build();

        // let header_unmapped = sam::Header::builder()
        //     .set_header(header_input.header().expect("No header present").clone())
        //     .add_reference_sequence("unmapped",  Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1e6 as usize)?)) // Provide a dummy reference sequence argument
        //     .build();

        let header_unmapped = sam::Header::builder()
            .set_header(header_input.header().expect("No header present").clone())
            .set_reference_sequences(reference_seqs)
            .build();

        Ok(BamHeaders {
            header_input,
            header_endogenous,
            header_exogenous,
            header_both_genomes,
            header_unmapped,
        })
    }

    fn write_headers(&mut self, headers: &BamHeaders) -> Result<()> {

        self.bam_endogenous
            .write_header(&headers.header_endogenous)?;
        self.bam_exogenous.write_header(&headers.header_exogenous)?;
        self.bam_both_genomes
            .write_header(&headers.header_both_genomes)?;
        self.bam_unmapped.write_header(&headers.header_unmapped)?;
        Ok(())
    }

    pub fn split(mut self, exogenous_prefix: &[u8]) -> Result<SplitStats> {
        let headers = self.make_headers(exogenous_prefix).expect("Error making headers");
        self.write_headers(&headers)?;

        let mut stats = SplitStats::new("SplitBam".to_string());

        let chromosomes: Vec<String> = headers
            .header_input
            .reference_sequences()
            .keys()
            .cloned()
            .map(|x| x.to_string())
            .collect();

        let queue_endogenous = ArrayQueue::new(1000);
        let queue_exogenous = ArrayQueue::new(1000);
        let queue_both_genomes = ArrayQueue::new(1000);
        let queue_unmapped = ArrayQueue::new(1000);


        // For each output channel, spawn a thread that writes records to the appropriate BAM file.
        let outputs = vec![
            (
                self.bam_endogenous,
                queue_endogenous,
                headers.header_endogenous.clone(),
            ),
            (
                self.bam_exogenous,
                queue_exogenous,
                headers.header_exogenous.clone(),
            ),
            (
                self.bam_both_genomes,
                queue_both_genomes,
                headers.header_both_genomes.clone(),
            ),
            (
                self.bam_unmapped,
                queue_unmapped,
                headers.header_unmapped.clone(),
            ),
        ];

        

        // Spawn threads to write records to the appropriate BAM file.
        for (writer, queue, header) in outputs {
            let _ = std::thread::spawn(move || {
                let mut writer = writer;
                writer.write_header(&header).expect("Error writing header");

                while let Some(rx) = queue.pop() {
                    writer.write_record(&header, &rx).expect("Error writing record");
                }
                
        });
        }
        
        
        let bam_input_path = self.bam_input.clone();

        // Spawn a thread for each chromosome to read records from the input BAM file and send them to the appropriate output channel.
        chromosomes.into_par_iter().for_each(|chromosome|
 {
                        
            let mut reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(&bam_input_path)
                .expect("Error building reader");
            
            let header = reader.read_header().expect("Error reading header");
            let mut stats = SplitStats::new(chromosome.clone());

            let chrom = chromosome.clone();
            let start = "1".to_string();
            let end  = header.reference_sequences().get(chrom.as_bytes()).expect("Error getting reference sequence").length().to_string();
            let coords = format!("{}:{}-{}", chrom, start, end);
            let region = &coords.parse().expect("Error parsing region");
            let query = reader.query(&header, region).expect("Failed to get query");

            for record in query {
                let record = record.expect("Error reading record");
                let record_clone = record.clone();

                if record.flags().is_unmapped() {
                    unmapped_tx
                        .try_send(record).unwrap();
                    stats.add_unmapped();
                    continue;
                } else if record.flags().is_qc_fail() {
                    unmapped_tx
                        .try_send(record).unwrap();
                    stats.add_qcfail();
                    continue;
                } else if record.flags().is_duplicate() {
                    unmapped_tx
                        .try_send(record).unwrap();
                    stats.add_duplicate();
                    continue;
                } else if record.flags().is_secondary() {
                    unmapped_tx
                        .try_send(record).unwrap();
                    stats.add_secondary();
                    continue;
                } else if record.mapping_quality().expect("No mapping quality").get() < 30 {
                    unmapped_tx
                        .try_send(record).unwrap();
                    stats.add_low_maq();
                    continue;
                } else if !record.flags().is_mate_unmapped() {
                    let r1_seq_id = record
                        .reference_sequence_id()
                        .expect("No reference sequence ID")
                        .expect("Failed to get reference sequence ID");
                    let r1_seq_name = headers
                        .header_input
                        .reference_sequences()
                        .get_index(r1_seq_id)
                        .expect("Failed to get reference sequence name")
                        .0;
                    let r2_seq_id = record
                        .mate_reference_sequence_id()
                        .expect("No mate reference sequence ID")
                        .expect("Failed to get mate reference sequence ID");
                    let r2_seq_name = headers
                        .header_input
                        .reference_sequences()
                        .get_index(r2_seq_id)
                        .expect("Failed to get mate reference sequence name")
                        .0;

                    if r1_seq_name.starts_with(exogenous_prefix)
                        && r2_seq_name.starts_with(exogenous_prefix)
                    {
                        exogenous_tx
                            .try_send(record).unwrap();
                        stats.add_exogenous();
                        continue;
                    } else if r1_seq_name.starts_with(exogenous_prefix)
                        || r2_seq_name.starts_with(exogenous_prefix)
                    {
                        both_genomes_tx
                            .try_send(record).unwrap();
                        stats.add_both_genomes();
                        continue;
                    } else {
                        endogenous_tx
                            .try_send(record).unwrap();
                    }
                }
            }
        });


        Ok(stats)
    }
}
