use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{error, info};
use std::process::exit;
use std::{
    io::Write,
    path::{Path, PathBuf},
};

mod count;
mod filter;
mod intervals;
mod spikein;
mod utils;
// mod meta;

use crate::utils::FileType;


pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}


#[derive(Parser)]
#[command(author, version, about, long_about = None, styles=get_styles())]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    /// Verbosity level
    #[arg(short, long, required = false, default_value = "0")]
    verbose: u8,

    // Filter parameters
    /// Properly paired reads only
    #[arg(long, action = clap::ArgAction::SetTrue)]
    proper_pair: bool,

    /// Minimum mapping quality
    #[arg(long, required = false, default_value = "20")]
    min_mapq: u8,

    /// Minimum read length
    #[arg(long, required = false, default_value = "20")]
    min_length: u32,

    /// Maximum read length
    #[arg(long, required = false, default_value = "1000")]
    max_length: u32,

    /// Blacklisted locations in BED format
    #[arg(long, required = false)]
    blacklisted_locations: Option<PathBuf>,

    /// Whitelisted barcodes in a text file (one barcode per line)
    #[arg(long, required = false)]
    whitelisted_barcodes: Option<PathBuf>,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate coverage from a BAM file and write to a bedGraph or bigWig file
    BamCoverage {
        /// Bam file for processing
        #[arg(short, long)]
        bam: PathBuf,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Bin size for coverage calculation
        #[arg(short = 's', long)]
        bin_size: Option<u64>,

        /// Normalization method to use
        #[arg(short, long, default_value = "raw")]
        norm_method: Option<utils::NormalizationMethod>,

        /// Scaling factor for the pileup
        #[arg(short = 'f', long)]
        scale_factor: Option<f32>,

        /// Use the fragment or the read for counting
        #[arg(long, action = clap::ArgAction::SetTrue)]
        use_fragment: bool,
    },

    /// Split a BAM file into endogenous and exogenous reads
    SplitExogenous {
        /// Input BAM file
        #[arg(short, long)]
        input: PathBuf,

        /// Output prefix
        #[arg(short, long)]
        output: PathBuf,

        /// Prefix for exogenous sequences
        #[arg(short, long)]
        exogenous_prefix: String,

        /// Path for stats output
        #[arg(short, long)]
        stats: Option<PathBuf>,
    },
}

fn main() {
    colog::init();

    let cli = Cli::parse();

    let whitelisted_barcodes = match &cli.whitelisted_barcodes {
        Some(whitelist) => {
            let barcodes = utils::CellBarcodes::from_csv(whitelist).expect("Failed to read barcodes");
            Some(barcodes.barcodes())
        }
        None => None,
    };

    match &cli.command {
        Some(Commands::BamCoverage {
            bam,
            output,
            bin_size,
            norm_method,
            scale_factor,
            use_fragment,
        }) => {
            if !bam.exists() {
                error!("BAM file does not exist");
                exit(1);
            } else if !bam.with_extension("bam.bai").exists() {
                error!("BAM index file does not exist. Please create the index file using samtools index command");
            }



            let filter = filter::BamReadFilter::new(
                cli.proper_pair,
                Some(cli.min_mapq),
                Some(cli.min_length),
                Some(cli.max_length),
                None,
                whitelisted_barcodes,
            );
            let coverage = count::BamPileup::new(
                bam.clone(),
                bin_size.unwrap_or(50),
                norm_method
                    .clone()
                    .unwrap_or(utils::NormalizationMethod::Raw),
                scale_factor.unwrap_or(1.0),
                *use_fragment,
                filter,
            );

            let outfile = match output {
                Some(output) => output.clone(),
                None => {
                    let mut output = bam.clone();
                    output.set_extension("bedgraph");
                    output
                }
            };

            let output_type = match outfile.extension() {
                Some(ext) => match ext.to_str() {
                    Some(ext) => match FileType::from_str(ext) {
                        Ok(file_type) => Some(file_type),
                        Err(e) => {
                            error!("{}", e);
                            None
                        }
                    },
                    None => {
                        error!("Could not determine file type from extension");
                        None
                    }
                },
                None => {
                    error!("Could not determine file type from extension");
                    None
                }
            };

            let result = match output_type {
                Some(FileType::Bedgraph) => coverage.to_bedgraph(outfile),
                Some(FileType::Bigwig) => {
                    // Check that bedGraphToBigWig is in the PATH

                    // Check if the bedGraphToBigWig tool is in the PATH
                    let bdg_to_bw = std::process::Command::new("bedGraphToBigWig").output();
                    match bdg_to_bw {
                        std::result::Result::Ok(_) => {
                            log::info!("bedGraphToBigWig tool found in PATH");
                        }
                        Err(_) => {
                            log::error!("bedGraphToBigWig tool not found in PATH. Please install the tool and try again");
                            std::process::exit(1);
                        }
                    };

                    coverage.to_bigwig(outfile)
                }
                None => {
                    exit(1);
                }
            };

            match result {
                Ok(_) => {
                    info!("Successfully wrote output");
                }
                Err(e) => {
                    error!("Error writing output: {}", e);
                    exit(1);
                }
            }
        }
        Some(Commands::SplitExogenous {
            input,
            output,
            exogenous_prefix,
            stats,
        }) => {
            let mut split = spikein::BamSplitter::new(
                input.to_path_buf(),
                output.to_path_buf(),
                exogenous_prefix.to_string(),
            )
            .expect("Failed to create BamSplitter");
            match split.split() {
                Ok(_) => {
                    info!("Successfully split BAM file");

                    if stats.is_none() {
                        return;
                    } else {
                        let stats = stats.as_ref().unwrap();
                        let split_stats = split.stats();
                        let json =
                            serde_json::to_string(&split_stats).expect("Failed to serialize stats");
                        let mut stats_file =
                            std::fs::File::create(stats).expect("Failed to create stats file");
                        stats_file
                            .write_all(json.as_bytes())
                            .expect("Failed to write stats file");
                        info!("Successfully wrote stats file to {}", stats.display());
                    }
                }
                Err(e) => {
                    error!("Error splitting BAM file: {}", e);
                    exit(1);
                }
            }
        }
        None => {
            eprintln!("No subcommand provided");
            std::process::exit(1);
        }
    }
}
