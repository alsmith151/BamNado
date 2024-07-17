use polars::prelude::*;
use ahash::AHashSet as HashSet;
use anyhow::Result;
use std::path::Path;

pub const CB: [u8; 2] = [b'C', b'B'];


pub struct Barcodes {
    barcodes: HashSet<String>,
}

impl Barcodes {
    pub fn new() -> Self {
        Self {
            barcodes: HashSet::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.barcodes.is_empty()
    }

    pub fn from_csv(file_path: &str) -> Result<Self> {

        let path = Path::new(file_path).to_path_buf();

        let df = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path))?
            .finish()?;
        let mut barcodes = HashSet::new();
        
        for barcode in df.column("barcode").unwrap().str().unwrap() {
            let barcode = barcode.unwrap().to_string();
            barcodes.insert(barcode);
        }

        println!("Number of barcodes: {}", barcodes.len());
        println!("First 10 barcodes: {:?}", barcodes.iter().take(10).collect::<Vec<&String>>());

        Ok(Self { barcodes })
    }

    pub fn len(&self) -> usize {
        self.barcodes.len()
    }

    pub fn barcodes_list(&self) -> Vec<String> {
        self.barcodes.iter().map(|x| x.to_string()).collect()
    }

    pub fn barcodes(&self) -> HashSet<String>{
        self.barcodes.clone()
    }
}



pub fn progress_bar(length: u64, message: String) -> indicatif::ProgressBar {
    // Progress bar
    let progress_style = indicatif::ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta}) {msg}",
    )
    .unwrap()
    .progress_chars("##-");

    let progress_bar = indicatif::ProgressBar::new(length as u64);
    progress_bar.set_style(progress_style.clone());
    progress_bar.set_message(message);
    progress_bar
}