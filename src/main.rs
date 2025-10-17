use anyhow::{Context, Result};
use bamnado::bam_utils::FileType;
use clap::{Parser, Subcommand};
use itertools::Itertools;
use log::info;
use noodles::bam;
use noodles::core::{Position, Region};
use rayon::prelude::*;
use std::str::FromStr;
use std::{
    io::Write,
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

    Fragments {
        /// Input BAM file
        #[arg(short, long)]
        input: PathBuf,

        /// Output prefix
        #[arg(short, long)]
        output: PathBuf,

        /// Minimum fragment length
        #[arg(long, default_value = "50")]
        min_fragment_length: u32,

        /// Maximum fragment length
        #[arg(long, default_value = "1000")]
        max_fragment_length: u32,

        /// Number of reads to process before collapsing and writing output to file toreduce memory usage
        #[arg(long, default_value = "1000000")]
        chunk_size: usize,
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
            let lapper = bamnado::bam_utils::lapper_chrom_name_to_lapper_chrom_id(lapper, stats)
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
                        let lapper = bamnado::bam_utils::lapper_chrom_name_to_lapper_chrom_id(
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

        Commands::Fragments {
            input,
            output,
            min_fragment_length,
            max_fragment_length,
            chunk_size,
        } => {

            // Validate input BAM file
            validate_bam_file(input)?;

            let header = bamnado::bam_utils::get_bam_header(input)?;
            let bam_stats = bamnado::BamStats::new(input.to_path_buf())
                .context("Failed to create BAM stats")?;
            let chromsizes = bam_stats.chromsizes_ref_name()?;

            // Create fragment structs
            let mut config = bamnado::fragments::FragmentConfig::default();
            config.min_distance = *min_fragment_length as i32;
            config.max_distance = *max_fragment_length as i32;

            info!("Fragment extraction configuration: {:?}", config);
            // Create fragment processor

            // Process BAM in chunks per chromosome
            chromsizes.clone().into_par_iter().try_for_each(
                |(chromosome, size)| -> anyhow::Result<()> {

                    //Create output directory if it does not exist
                    if let Some(parent) = output.parent() {
                        std::fs::create_dir_all(parent).context("Failed to create output directory")?;
                    }

                    // Create outfile path
                    let outpath =
                                format!("{}_{}.fragments.tsv.gz", output.display(), chromosome);
                    let outpath = PathBuf::from(outpath);


                    // Create a fragment tracker
                    let mut tracker = bamnado::fragments::FragmentTracker::new(config.clone());
                    let collapser = bamnado::fragments::FragmentCollapser::new(config.clone());

                    // Create a new reader for each thread
                    let mut reader = bam::io::indexed_reader::Builder::default()
                        .build_from_path(input)
                        .context(format!("Failed to open BAM file {}", input.display()))?;

                    let start = Position::try_from(1)?;
                    let end = Position::try_from(size as usize)?;

                    let records = reader
                        .query(&header, &Region::new(chromosome.clone(), start..=end))
                        .context("Failed to query BAM records")?;


                    records.chunks(*chunk_size).into_iter().try_for_each(
                        |chunk| -> anyhow::Result<()> {
                            // collect records in this chunk
                            let records: Result<Vec<_>, _> = chunk.collect();
                            let records = records.context("Failed to read BAM records")?;

                            // ensure there is at least one record in the chunk
                            let last = records
                                .last()
                                .ok_or_else(|| anyhow::anyhow!("No records in chunk"))?;

                            for (_, record) in records.iter().enumerate() {
                                tracker.process_record(record, &chromosome)?;
                            }

                            // compute current position based on the last record
                            let current_position: i32 = match last.alignment_start() {
                                Some(Ok(pos)) => {
                                    if last.flags().is_reverse_complemented() {
                                        // for reverse strand reads we use the alignment start
                                        usize::from(pos) as i32
                                    } else {
                                        // for forward strand reads use alignment start + read length
                                        // be careful with usize -> i32 conversion
                                        let read_len = last.sequence().len() as i32;
                                        usize::from(pos) as i32 + read_len
                                    }
                                }
                                Some(Err(e)) => {
                                    return Err(anyhow::anyhow!(
                                        "Failed to get alignment position: {}",
                                        e
                                    ));
                                }
                                None => {
                                    return Err(anyhow::anyhow!(
                                        "Alignment position is missing for a record"
                                    ));
                                }
                            };


                            // harvest complete fragments up to the current position
                            let fragments_complete =
                                tracker.harvest_complete_fragments(current_position);

                            
                            let fragments_collapsed = collapser.collapse(fragments_complete);
                            
                            bamnado::fragments::write_fragments(
                                &fragments_collapsed,
                                &outpath,
                                true,
                            )?;


                            Ok(())
                        },
                    )?;


                    // After processing all chunks, harvest any remaining fragments
                    let remaining_fragments = tracker.harvest_all();

                    let remaining_collapsed = collapser.collapse(remaining_fragments);
                    bamnado::fragments::write_fragments(
                        &remaining_collapsed,
                        &outpath,
                        true,
                    )?;

                    // Log stats after processing all chunks
                    let stats = tracker.stats();
                    info!(
                        "Chromosome: {}, Total fragments: {}, Incomplete fragments: {}",
                        chromosome, stats.1, stats.2
                    );

                    Ok(())
                },
            )?;

            // Merge per-chromosome fragment files into a single output file
            let mut final_outfile = std::fs::File::create(output)
                .context("Failed to create final fragments output file")?;

            // Just read the contents of each per-chromosome file and write to final output
            chromsizes.into_iter().try_for_each(|(chromosome, _size)| -> anyhow::Result<()> {
                let per_chrom_path =
                    format!("{}_{}.fragments.tsv.gz", output.display(), chromosome);
                let per_chrom_path = PathBuf::from(per_chrom_path);

                let mut per_chrom_file = std::fs::File::open(&per_chrom_path)
                    .context("Failed to open per-chromosome fragments file")?;

                std::io::copy(&mut per_chrom_file, &mut final_outfile)
                    .context("Failed to write to final fragments output file")?;

                // Remove the per-chromosome file
                std::fs::remove_file(&per_chrom_path)
                    .context("Failed to remove temporary per-chromosome fragments file")?;

                Ok(())
            }).expect("Encountered an error when making fragments");
            



        }
    }

    Ok(())
}

