use anyhow::{Context, Result};
use bamnado::bam_utils::FileType;
use clap::{Parser, Subcommand};
use log::info;
use polars::prelude::*;
use std::str::FromStr;
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

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

// Define common filter options as a trait to avoid code duplication
#[derive(Parser, Clone)]
struct FilterOptions {
    /// Filter reads based on strand
    #[arg(long, default_value = "both")]
    strand: bamnado::Strand,

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

    /// Selected read group
    #[arg(long, required = false)]
    read_group: Option<String>,

    /// Filter reads by SAM tag (e.g., VP, AS, etc.)
    #[arg(long, required = false)]
    filter_tag: Option<String>,

    /// Value for the filter tag (e.g., BCL2, TP, etc.)
    #[arg(long, required = false)]
    filter_tag_value: Option<String>,
}

// Common coverage options for both single and multi-BAM operations
#[derive(Parser, Clone)]
struct CoverageOptions {
    /// Bin size for coverage calculation
    #[arg(short = 's', long)]
    bin_size: Option<u64>,

    /// Normalization method to use
    #[arg(short, long, default_value = "raw")]
    norm_method: Option<bamnado::signal_normalization::NormalizationMethod>,

    /// Scaling factor for the pileup
    #[arg(short = 'f', long)]
    scale_factor: Option<f32>,

    /// Use the fragment or the read for counting
    #[arg(long, action = clap::ArgAction::SetTrue)]
    use_fragment: bool,

    /// Shift options for the pileup
    #[arg(long, default_value = "0,0,0,0")]
    shift: Option<bamnado::genomic_intervals::Shift>,

    /// Truncate options for the pileup
    #[arg(long)]
    truncate: Option<bamnado::genomic_intervals::Truncate>,

    /// Ignore scaffold chromosomes
    #[arg(long, action = clap::ArgAction::SetTrue)]
    ignore_scaffold: bool,
}

#[derive(Parser)]
#[command(author, version, about, long_about = None, styles=get_styles())]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Verbosity level
    #[arg(short, long, required = false, default_value = "2")]
    verbose: u8,
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

        #[command(flatten)]
        coverage_options: CoverageOptions,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Calculate coverage from multiple BAM files and write to a bedGraph or bigWig file
    MultiBamCoverage {
        /// List of BAM files for processing
        #[arg(short, long)]
        bams: Vec<PathBuf>,

        /// Output file name
        #[arg(short, long)]
        output: Option<PathBuf>,

        #[command(flatten)]
        coverage_options: CoverageOptions,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Compare two BigWig files
    #[command(name = "bigwig-compare")]
    CompareBigWigs {
        /// Path to the first BigWig file
        #[arg(long)]
        bw1: PathBuf,

        /// Path to the second BigWig file
        #[arg(long)]
        bw2: PathBuf,

        /// Output BigWig file path
        #[arg(short, long)]
        output: PathBuf,

        /// Comparison method
        #[arg(short, long, value_enum)]
        comparison: bamnado::bigwig_compare::Comparison,

        /// Bin size for comparison
        #[arg(short = 's', long, default_value = "50")]
        bin_size: u32,

        /// Chunk size for processing
        #[arg(long)]
        chunk_size: Option<usize>,

        #[arg(long)]
        pseudocount: Option<f64>,
    },

    /// Aggregate multiple BigWig files into one
    #[command(name = "bigwig-aggregate")]
    AggregateBigWigs {
        /// Paths to the BigWig files to aggregate
        #[arg(long, num_args = 1..)]
        bigwigs: Vec<PathBuf>,

        /// Output BigWig file path
        #[arg(short, long)]
        output: PathBuf,

        /// Aggregation method
        #[arg(short, long, value_enum)]
        method: bamnado::bigwig_compare::AggregationMode,

        /// Bin size for aggregation
        #[arg(short = 's', long, default_value = "50")]
        bin_size: u32,

        /// Chunk size for processing
        #[arg(long)]
        chunk_size: Option<usize>,

        /// Pseudocount value to add to all values
        #[arg(long)]
        pseudocount: Option<f64>,
    },

    /// Collapse adjacent equal-score bins in a BedGraph
    #[command(name = "collapse-bedgraph")]
    CollapseBedgraph {
        /// Input BedGraph file (defaults to stdin)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output BedGraph file (defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Split a BAM file based on a set of defined filters
    Split {
        /// Input BAM file
        #[arg(short, long)]
        input: PathBuf,

        /// Output prefix
        #[arg(short, long)]
        output: PathBuf,

        #[command(flatten)]
        filter_options: FilterOptions,
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

        /// Allow unknown MAPQ values - useful for BAM files with MAPQ=255 i.e. STAR generated BAM files
        #[arg(long, action = clap::ArgAction::SetTrue)]
        allow_unknown_mapq: bool,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Filter and/or adjust reads in a BAM file e.g. shift for Tn5 insertion
    Modify {
        /// Input BAM file
        #[arg(short, long)]
        input: PathBuf,

        /// Output prefix
        #[arg(short, long)]
        output: PathBuf,

        /// Reads to ignore during modification
        #[command(flatten)]
        filter_options: FilterOptions,

        #[arg(long, action = clap::ArgAction::SetTrue)]
        tn5_shift: bool,
    },
}

// Helper functions to reduce code duplication
fn validate_bam_file(bam_path: &Path) -> Result<()> {
    if !bam_path.exists() {
        return Err(anyhow::anyhow!(
            "BAM file does not exist: {}",
            bam_path.display()
        ));
    }

    let index_path = bam_path.with_extension("bam.bai");
    if !index_path.exists() {
        return Err(anyhow::anyhow!(
            "BAM index file does not exist for {}. Please create the index file using samtools index command",
            bam_path.display()
        ));
    }

    Ok(())
}

fn create_filter_from_options(
    filter_options: &FilterOptions,
    bam_stats: Option<&bamnado::BamStats>,
) -> Result<bamnado::read_filter::BamReadFilter> {
    // Process whitelisted barcodes
    let whitelisted_barcodes = match &filter_options.whitelisted_barcodes {
        Some(whitelist) => {
            let barcodes =
                bamnado::CellBarcodes::from_csv(whitelist).context("Failed to read barcodes")?;
            Some(barcodes.barcodes())
        }
        None => None,
    };

    // Process blacklisted locations
    let blacklisted_locations = match (&filter_options.blacklisted_locations, bam_stats) {
        (Some(blacklist), Some(stats)) => {
            let lapper = bamnado::bam_utils::bed_to_lapper(blacklist.clone())
                .context("Failed to read blacklisted locations")?;

            let lapper = bamnado::bam_utils::convert_lapper_chrom_names_to_ids(lapper, stats)
                .context("Failed to convert chrom names to chrom ids")?;

            // Merge overlapping intervals
            let lapper = lapper
                .into_iter()
                .map(|(chrom, mut intervals)| {
                    intervals.merge_overlaps();
                    (chrom, intervals)
                })
                .collect();

            Some(lapper)
        }
        (Some(_blacklist), None) => {
            return Err(anyhow::anyhow!(
                "Blacklisted locations provided but no BAM stats available"
            ));
        }
        _ => None,
    };

    println!("Blacklisted locations: {blacklisted_locations:?}");

    Ok(bamnado::read_filter::BamReadFilter::new(
        filter_options.strand.into(),
        filter_options.proper_pair,
        Some(filter_options.min_mapq),
        Some(filter_options.min_length),
        Some(filter_options.max_length),
        filter_options.read_group.clone(),
        blacklisted_locations,
        whitelisted_barcodes,
        filter_options.filter_tag.clone(),
        filter_options.filter_tag_value.clone(),
    ))
}

fn process_output_file_type(output: &Path) -> Result<FileType> {
    match output.extension() {
        Some(ext) => match ext.to_str() {
            Some(ext) => FileType::from_str(ext)
                .map_err(anyhow::Error::msg)
                .context(format!("Unsupported file extension: {ext}")),
            None => Err(anyhow::anyhow!(
                "Could not determine file type from extension"
            )),
        },
        None => Err(anyhow::anyhow!("No file extension found")),
    }
}

fn main() -> Result<()> {
    colog::init();

    let cli = Cli::parse();
    let log_level = match cli.verbose {
        0 => log::LevelFilter::Error,
        1 => log::LevelFilter::Warn,
        2 => log::LevelFilter::Info,
        3 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    log::set_max_level(log_level);

    match &cli.command {
        Commands::BamCoverage {
            bam,
            output,
            coverage_options,
            filter_options,
        } => {
            // Validate BAM file
            validate_bam_file(bam)?;

            // Get the BAM stats
            let bam_stats =
                bamnado::BamStats::new(bam.clone()).context("Failed to read BAM stats")?;

            // Create filter
            let filter = create_filter_from_options(filter_options, Some(&bam_stats))?;

            // Create pileup
            let coverage = bamnado::coverage_analysis::BamPileup::new(
                bam.clone(),
                coverage_options.bin_size.unwrap_or(50),
                coverage_options
                    .norm_method
                    .clone()
                    .unwrap_or(bamnado::signal_normalization::NormalizationMethod::Raw),
                coverage_options.scale_factor.unwrap_or(1.0),
                coverage_options.use_fragment,
                filter,
                true,
                coverage_options.ignore_scaffold,
                coverage_options.shift,
                coverage_options.truncate,
            );

            // Determine output file
            let outfile = match output {
                Some(output) => output.clone(),
                None => {
                    let mut output = bam.clone();
                    output.set_extension("bedgraph");
                    output
                }
            };

            // Process output based on file type
            let output_type = process_output_file_type(&outfile)?;

            match output_type {
                FileType::Bedgraph => {
                    coverage
                        .to_bedgraph(outfile)
                        .context("Failed to write bedgraph output")?;
                }
                FileType::Bigwig => {
                    coverage
                        .to_bigwig(outfile)
                        .context("Failed to write bigwig output")?;
                }
                FileType::TSV => {
                    return Err(anyhow::anyhow!(
                        "TSV output is not supported for single BAM coverage"
                    ));
                }
            }

            info!("Successfully wrote output");
        }

        Commands::MultiBamCoverage {
            bams,
            output,
            coverage_options,
            filter_options,
        } => {
            // Validate all BAM files
            for bam in bams {
                validate_bam_file(bam)?;
            }

            // Determine output file
            let output = output.clone().unwrap_or_else(|| {
                let mut output = PathBuf::new();
                output.push("output");
                output.set_extension("tsv");
                output
            });

            // Get whitelisted barcodes for all BAMs
            let whitelisted_barcodes = match &filter_options.whitelisted_barcodes {
                Some(whitelist) => {
                    let barcodes = bamnado::bam_utils::CellBarcodesMulti::from_csv(whitelist)
                        .context("Failed to read barcodes")?;
                    Some(barcodes.barcodes())
                }
                None => None,
            };

            // Create filters and pileups for each BAM file
            let mut pileups = Vec::new();
            for (index, bam) in bams.iter().enumerate() {
                // Get BAM stats
                let bam_stats =
                    bamnado::BamStats::new(bam.clone()).context("Failed to read BAM stats")?;

                // Process blacklisted locations
                let blacklisted_locations = match &filter_options.blacklisted_locations {
                    Some(blacklist) => {
                        let lapper = bamnado::bam_utils::bed_to_lapper(blacklist.clone())
                            .context("Failed to read blacklisted locations")?;
                        let lapper = bamnado::bam_utils::convert_lapper_chrom_names_to_ids(
                            lapper, &bam_stats,
                        )
                        .context("Failed to convert chrom names to chrom ids")?;
                        Some(lapper)
                    }
                    None => None,
                };

                // Get specific whitelisted barcodes for this BAM
                let bam_barcodes = whitelisted_barcodes
                    .as_ref()
                    .map(|whitelist| whitelist[index].clone());

                // Create filter
                let filter = bamnado::read_filter::BamReadFilter::new(
                    filter_options.strand.into(),
                    filter_options.proper_pair,
                    Some(filter_options.min_mapq),
                    Some(filter_options.min_length),
                    Some(filter_options.max_length),
                    filter_options.read_group.clone(),
                    blacklisted_locations,
                    bam_barcodes,
                    filter_options.filter_tag.clone(),
                    filter_options.filter_tag_value.clone(),
                );

                // Create pileup
                pileups.push(bamnado::coverage_analysis::BamPileup::new(
                    bam.clone(),
                    coverage_options.bin_size.unwrap_or(50),
                    coverage_options
                        .norm_method
                        .clone()
                        .unwrap_or(bamnado::signal_normalization::NormalizationMethod::Raw),
                    coverage_options.scale_factor.unwrap_or(1.0),
                    coverage_options.use_fragment,
                    filter,
                    false,
                    coverage_options.ignore_scaffold,
                    coverage_options.shift,
                    coverage_options.truncate,
                ));
            }

            // Create multi-BAM pileup
            let coverage = bamnado::coverage_analysis::MultiBamPileup::new(pileups);

            // Process output based on file type
            let output_type = process_output_file_type(&output)?;

            match output_type {
                FileType::TSV => {
                    coverage
                        .to_tsv(&output)
                        .context("Failed to write TSV output")?;
                }
                _ => {
                    return Err(anyhow::anyhow!(
                        "Unsupported output format. Currently only TSV is supported for multi-BAM coverage"
                    ));
                }
            }

            info!("Successfully wrote output");
        }

        Commands::SplitExogenous {
            input,
            output,
            exogenous_prefix,
            stats,
            allow_unknown_mapq: _,
            filter_options,
        } => {
            // Validate input BAM file
            validate_bam_file(input)?;

            // Create filter
            let _filter = create_filter_from_options(filter_options, None)?;

            // Create and run BAM splitter
            let mut split = bamnado::spike_in_analysis::BamSplitter::new(
                input.clone(),
                output.clone(),
                exogenous_prefix.clone(),
            )
            .context("Failed to create BamSplitter")?;

            split.split().context("Failed to split BAM file")?;

            info!("Successfully split BAM file");

            // Write stats if requested
            if let Some(stats_path) = stats {
                let split_stats = split.stats();
                let json =
                    serde_json::to_string(&split_stats).context("Failed to serialize stats")?;

                let mut stats_file =
                    std::fs::File::create(stats_path).context("Failed to create stats file")?;

                stats_file
                    .write_all(json.as_bytes())
                    .context("Failed to write stats file")?;

                info!("Successfully wrote stats file to {}", stats_path.display());
            }
        }

        Commands::Split {
            input,
            output,
            filter_options,
        } => {
            // Validate input BAM file
            validate_bam_file(input)?;

            // Create filter
            let filter = create_filter_from_options(filter_options, None)?;

            // Create BAM splitter
            let split =
                bamnado::bam_splitter::BamFilterer::new(input.clone(), output.clone(), filter);

            // Run splitter asynchronously
            let rt = tokio::runtime::Builder::new_current_thread()
                .enable_all()
                .build()
                .context("Failed to build Tokio runtime")?;

            let split_future = split.split_async();
            rt.block_on(split_future)
                .context("Failed to split BAM file")?;

            info!("Successfully split BAM file");
        }

        Commands::Modify {
            input,
            output,
            filter_options,
            tn5_shift,
        } => {
            // Validate input BAM file
            validate_bam_file(input)?;

            // Create filter
            let filter = create_filter_from_options(filter_options, None)?;

            // Create BAM modifier
            let modifier = bamnado::bam_modifier::BamModifier::new(
                input.clone(),
                output.clone(),
                filter,
                *tn5_shift,
            );

            // Run modifier asynchronously
            let rt = tokio::runtime::Builder::new_current_thread()
                .enable_all()
                .build()
                .context("Failed to build Tokio runtime")?;

            let modify_future = modifier.run();
            rt.block_on(modify_future)
                .context("Failed to modify BAM file")?;

            info!("Successfully modified BAM file");
        }

        Commands::CompareBigWigs {
            bw1,
            bw2,
            output,
            comparison,
            bin_size,
            chunk_size,
            pseudocount,
        } => {
            bamnado::bigwig_compare::compare_bigwigs(
                bw1,
                bw2,
                output,
                comparison.clone(),
                *bin_size,
                *chunk_size,
                *pseudocount,
            )
            .context("Failed to compare BigWig files")?;

            info!(
                "Successfully compared BigWig files and wrote output to {}",
                output.display()
            );
        }

        Commands::AggregateBigWigs {
            bigwigs,
            output,
            method,
            bin_size,
            chunk_size,
            pseudocount,
        } => {
            if bigwigs.is_empty() {
                return Err(anyhow::anyhow!("At least one BigWig file must be provided"));
            }

            bamnado::bigwig_compare::aggregate_bigwigs(
                bigwigs,
                output,
                method.clone(),
                *bin_size,
                *chunk_size,
                *pseudocount,
            )
            .context("Failed to aggregate BigWig files")?;

            info!(
                "Successfully aggregated {} BigWig files and wrote output to {}",
                bigwigs.len(),
                output.display()
            );
        }
        Commands::CollapseBedgraph { input, output } => {
            // read bedgraph
            let reader: Box<dyn BufRead> = match input {
                Some(p) => Box::new(BufReader::new(File::open(p)?)),
                None => Box::new(BufReader::new(io::stdin())),
            };

            let mut chroms: Vec<String> = Vec::new();
            let mut starts: Vec<i64> = Vec::new();
            let mut ends: Vec<i64> = Vec::new();
            let mut scores: Vec<f64> = Vec::new();

            for line in reader.lines() {
                let l = line?;
                let l = l.trim();
                if l.is_empty() || l.starts_with('#') {
                    continue;
                }
                let parts: Vec<&str> = l.split_whitespace().collect();
                if parts.len() < 4 {
                    continue;
                }
                chroms.push(parts[0].to_string());
                starts.push(parts[1].parse()?);
                ends.push(parts[2].parse()?);
                scores.push(parts[3].parse()?);
            }

            let df = DataFrame::new(vec![
                Column::new("chrom".into(), chroms),
                Column::new("start".into(), starts),
                Column::new("end".into(), ends),
                Column::new("score".into(), scores),
            ])?;

            let collapsed = bamnado::bedgraph_utils::collapse_adjacent_bins(df)?;

            match output {
                Some(p) => {
                    let f = File::create(p)?;
                    let mut out_df = collapsed;
                    let w = CsvWriter::new(f);
                    w.include_header(false)
                        .with_separator(b'\t')
                        .finish(&mut out_df)?;
                }
                None => {
                    let stdout = io::stdout();
                    let handle = stdout.lock();
                    let mut out_df = collapsed;
                    let w = CsvWriter::new(handle);
                    w.include_header(false)
                        .with_separator(b'\t')
                        .finish(&mut out_df)?;
                }
            }
        }
    }

    Ok(())
}
