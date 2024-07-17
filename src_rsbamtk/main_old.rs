use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::{PathBuf};

pub mod shift;
pub mod filter;
pub mod split_exogenous;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    Shift {
        /// Bam file for processing
        #[arg(short, long)]
        bam: Option<PathBuf>,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    Subtract {
        /// Bed file for processing
        #[arg(short='r', long="regions")]
        regions: Option<PathBuf>,

        /// Bam file for processing
        #[arg(short='b', long="bam")]
        bam: Option<PathBuf>,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Number of threads to use
        #[arg(short, long)]
        threads: Option<usize>,
    },

    Split {
        /// Bam file for processing
        #[arg(short, long)]
        bam: PathBuf,

        /// Prefix to use for exogenous spike-in reads
        /// If not provided will default to dm6_
        #[arg(short, long)]
        exogenous_prefix: Option<String>,


        /// Output file prefix. The output files will be named as prefix_X.bam
        #[arg(short, long)]
        output: Option<PathBuf>,

    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Shift { bam, output }) => match (bam, output) {
            (Some(bam_file), Some(output_file)) => {
                shift::atac_shift_bam(bam_file, output_file).with_context(|| {
                    format!(
                        "Shifting reads failed for file `{}`",
                        bam_file.to_string_lossy()
                    )
                })?;
            }
            (Some(bam_file), None) => {
                let output_file = &PathBuf::from("shifted.bam");
                shift::atac_shift_bam(bam_file, output_file).with_context(|| {
                    format!(
                        "Shifting reads failed for file `{}`",
                        bam_file.to_string_lossy()
                    )
                })?
            }
            _ => {
                println!("Options not provided, will not run")
            }
        },

        Some(Commands::Subtract {
            regions: bed,
            bam,
            output,
            threads,
        }) => {
            println!("Running subtract subcommand. Will subtract regions from BAM file.");

            match (bed, bam) {
                (Some(bed_file), Some(bam_file)) => {
                    let output = match output {
                        Some(output) => output.to_owned(),
                        None => PathBuf::from("subtracted.bam"),
                    };
                    let threads = match threads {
                        Some(threads) => *threads,
                        None => 1,
                    };

                    println!("BED file: {}", bed_file.to_string_lossy());
                    println!("BAM file: {}", bam_file.to_string_lossy());
                    println!("Output file: {}", output.to_string_lossy());
                    println!("Threads: {}", threads);
                    
                    filter::remove_regions_from_bam(
                        bed_file.to_path_buf(),
                        bam_file.to_path_buf(),
                        output,
                        threads,
                    )?;
                }
                _ => {
                    println!("Options not provided, will not run");
                    return Ok(());
                }
            }
        }

        Some(Commands::Split { bam, exogenous_prefix, output }) => match (bam, output) {
            (bam_file, Some(output_file)) => {
                let exogenous_prefix = match exogenous_prefix {
                    Some(prefix) => prefix.to_owned(),
                    None => "dm6_".to_string(),
                };
                let mut  splitter =  split_exogenous::SplitBam::new(bam_file.to_path_buf(), output_file.to_path_buf())?;
                let stats = splitter.split(exogenous_prefix.as_bytes())?;

                stats.print();
            }

            _ => {
                println!("Options not provided, will not run")
            }
        },

        _ => {
            println!("Subcommand not provided, will not run")
        }
    }
    Ok(())
}
