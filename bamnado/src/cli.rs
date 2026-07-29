use crate::bam_utils::FileType;
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use comfy_table::{
    Cell, CellAlignment, ContentArrangement, Table, modifiers::UTF8_ROUND_CORNERS,
    presets::UTF8_FULL,
};
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
#[command(next_help_heading = "Read Filters")]
struct FilterOptions {
    /// Count only reads from the selected strand.
    #[arg(long, default_value = "both", value_name = "STRAND")]
    strand: crate::Strand,

    /// Keep only properly paired reads.
    #[arg(
        long = "proper-pairs",
        visible_alias = "proper-pair",
        action = clap::ArgAction::SetTrue
    )]
    proper_pair: bool,

    /// Minimum mapping quality.
    #[arg(long, default_value = "20", value_name = "INT")]
    min_mapq: u8,

    /// Minimum read length in base pairs.
    #[arg(long, default_value = "20", value_name = "BP")]
    min_length: u32,

    /// Maximum read length in base pairs.
    #[arg(long, default_value = "1000", value_name = "BP")]
    max_length: u32,

    /// BED file of regions to exclude.
    #[arg(
        long = "blacklist",
        visible_alias = "blacklisted-locations",
        value_name = "BED"
    )]
    blacklisted_locations: Option<PathBuf>,

    /// Text file of barcodes to keep, one barcode per line.
    #[arg(
        long = "barcode-allowlist",
        visible_alias = "whitelisted-barcodes",
        value_name = "FILE"
    )]
    whitelisted_barcodes: Option<PathBuf>,

    /// Keep only reads from this read group.
    #[arg(long, value_name = "RG")]
    read_group: Option<String>,

    /// Keep reads with this SAM tag.
    #[arg(long = "tag", visible_alias = "filter-tag", value_name = "TAG")]
    filter_tag: Option<String>,

    /// Required value for `--tag`.
    #[arg(
        long = "tag-value",
        visible_alias = "filter-tag-value",
        value_name = "VALUE"
    )]
    filter_tag_value: Option<String>,

    /// Minimum fragment length (TLEN) in base pairs.
    #[arg(
        long = "min-fragment-len",
        visible_alias = "min-fragment-length",
        value_name = "BP"
    )]
    min_fragment_length: Option<u32>,

    /// Maximum fragment length (TLEN) in base pairs.
    #[arg(
        long = "max-fragment-len",
        visible_alias = "max-fragment-length",
        value_name = "BP"
    )]
    max_fragment_length: Option<u32>,

    /// Exclude PCR/optical duplicate reads.
    #[arg(
        long = "ignore-duplicates",
        visible_alias = "no-duplicates",
        action = clap::ArgAction::SetTrue
    )]
    ignore_duplicates: bool,

    /// Exclude secondary alignments.
    #[arg(long = "ignore-secondary", action = clap::ArgAction::SetTrue)]
    ignore_secondary: bool,

    /// Exclude supplementary alignments.
    #[arg(long = "ignore-supplementary", action = clap::ArgAction::SetTrue)]
    ignore_supplementary: bool,
}

// Common coverage options for both single and multi-BAM operations
#[derive(Parser, Clone)]
#[command(next_help_heading = "Coverage Options")]
struct CoverageOptions {
    /// Bin size for the output track.
    #[arg(short = 's', long, value_name = "BP")]
    bin_size: Option<u64>,

    /// Signal normalization method.
    #[arg(
        short = 'n',
        long = "normalize",
        visible_alias = "norm-method",
        default_value = "raw",
        value_name = "METHOD"
    )]
    norm_method: Option<crate::signal_normalization::NormalizationMethod>,

    /// Multiply the final signal by this factor.
    #[arg(short = 'f', long, value_name = "FLOAT")]
    scale_factor: Option<f32>,

    /// Count fragments instead of individual read alignments.
    #[arg(
        long = "fragment-counts",
        visible_alias = "use-fragment",
        action = clap::ArgAction::SetTrue
    )]
    use_fragment: bool,

    /// Shift read or fragment ends before pileup.
    #[arg(long, default_value = "0,0,0,0", value_name = "L,R,FL,FR")]
    shift: Option<crate::genomic_intervals::Shift>,

    /// Trim read or fragment ends before pileup.
    #[arg(long, value_name = "L,R,FL,FR")]
    truncate: Option<crate::genomic_intervals::Truncate>,

    /// Skip scaffold or unplaced chromosomes.
    #[arg(
        long = "ignore-scaffolds",
        visible_alias = "ignore-scaffold",
        action = clap::ArgAction::SetTrue
    )]
    ignore_scaffold: bool,

    /// Threads to use while writing BigWig output.
    #[arg(long, default_value = "6", value_name = "N")]
    threads: u32,
}

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Fast BAM and BigWig utilities for genomics workflows.",
    long_about = "Fast BAM and BigWig utilities for genomics workflows.\n\nUse `bamnado <command> --help` to see command-specific options and examples.",
    styles = get_styles(),
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Log verbosity from 0 (quiet) to 4 (debug).
    #[arg(short, long, default_value = "2", global = true, value_name = "LEVEL")]
    verbose: u8,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate a coverage track from one BAM file.
    #[command(
        name = "bam-coverage",
        visible_alias = "coverage",
        after_help = "Example:\n  bamnado bam-coverage --bam sample.bam --output sample.bw --bin-size 50 --normalize cpm"
    )]
    BamCoverage {
        /// Input BAM file.
        #[arg(short, long, value_name = "BAM")]
        bam: PathBuf,

        /// Output bedGraph or BigWig path.
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        #[command(flatten)]
        coverage_options: CoverageOptions,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Merge coverage from multiple BAM files into one track.
    #[command(name = "multi-bam-coverage", visible_alias = "multi-coverage")]
    MultiBamCoverage {
        /// Input BAM files.
        #[arg(short, long, value_name = "BAM")]
        bams: Vec<PathBuf>,

        /// Output bedGraph or BigWig path.
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        #[command(flatten)]
        coverage_options: CoverageOptions,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Compare two BigWig files bin by bin.
    #[command(name = "bigwig-compare", visible_alias = "compare-bigwigs")]
    CompareBigWigs {
        /// First BigWig input.
        #[arg(long, value_name = "BW")]
        bw1: PathBuf,

        /// Second BigWig input.
        #[arg(long, value_name = "BW")]
        bw2: PathBuf,

        /// Output BigWig path.
        #[arg(short, long, value_name = "BW")]
        output: PathBuf,

        /// Comparison method.
        #[arg(short, long, value_enum, value_name = "METHOD")]
        comparison: crate::bigwig_compare::Comparison,

        /// Bin size for comparison.
        #[arg(short = 's', long, default_value = "50", value_name = "BP")]
        bin_size: u32,

        /// Value added to both tracks before comparison.
        #[arg(long, value_name = "FLOAT")]
        pseudocount: Option<f64>,

        /// Scale factor applied to bw1 before comparison.
        #[arg(long, value_name = "FLOAT")]
        scale_factor_bw1: Option<f64>,

        /// Scale factor applied to bw2 before comparison.
        #[arg(long, value_name = "FLOAT")]
        scale_factor_bw2: Option<f64>,

        /// Number of threads to use for BigWig writing.
        #[arg(long, default_value = "6")]
        threads: u32,
    },

    /// Aggregate multiple BigWig files into one track.
    #[command(name = "bigwig-aggregate", visible_alias = "aggregate-bigwigs")]
    AggregateBigWigs {
        /// BigWig files to aggregate.
        #[arg(long, num_args = 1.., value_name = "BW")]
        bigwigs: Vec<PathBuf>,

        /// Output BigWig path.
        #[arg(short, long, value_name = "BW")]
        output: PathBuf,

        /// Aggregation method.
        #[arg(short, long, value_enum, value_name = "METHOD")]
        method: crate::bigwig_compare::AggregationMode,

        /// Bin size for aggregation.
        #[arg(short = 's', long, default_value = "50", value_name = "BP")]
        bin_size: u32,

        /// Value added to all inputs before aggregation.
        #[arg(long, value_name = "FLOAT")]
        pseudocount: Option<f64>,

        /// Per-file scale factors (one per --bigwigs input, in order).
        #[arg(long, num_args = 0.., value_name = "FLOAT")]
        scale_factors: Vec<f64>,

        /// Number of threads to use for BigWig writing.
        #[arg(long, default_value = "6")]
        threads: u32,
    },

    /// Collapse adjacent equal-score bins in a bedGraph.
    #[command(name = "collapse-bedgraph", visible_alias = "collapse")]
    CollapseBedgraph {
        /// Input bedGraph path. Reads from stdin if omitted.
        #[arg(short, long, value_name = "BEDGRAPH")]
        input: Option<PathBuf>,

        /// Output bedGraph path. Writes to stdout if omitted.
        #[arg(short, long, value_name = "BEDGRAPH")]
        output: Option<PathBuf>,
    },

    /// Split a BAM file using the supplied filters.
    Split {
        /// Input BAM file.
        #[arg(short, long, value_name = "BAM")]
        input: PathBuf,

        /// Output prefix.
        #[arg(short, long, value_name = "PREFIX")]
        output: PathBuf,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Split a BAM file into endogenous and exogenous reads.
    #[command(visible_alias = "split-spikein")]
    SplitExogenous {
        /// Input BAM file.
        #[arg(short, long, value_name = "BAM")]
        input: PathBuf,

        /// Output prefix.
        #[arg(short, long, value_name = "PREFIX")]
        output: PathBuf,

        /// Reference-name prefix used to identify exogenous sequences.
        #[arg(short, long, value_name = "PREFIX")]
        exogenous_prefix: String,

        /// Optional path for summary statistics output.
        #[arg(short, long, value_name = "FILE")]
        stats: Option<PathBuf>,

        #[command(flatten)]
        filter_options: FilterOptions,
    },

    /// Filter and/or adjust reads in a BAM file.
    Modify {
        /// Input BAM file.
        #[arg(short, long, value_name = "BAM")]
        input: PathBuf,

        /// Output prefix.
        #[arg(short, long, value_name = "PREFIX")]
        output: PathBuf,

        /// Reads to exclude before modification.
        #[command(flatten)]
        filter_options: FilterOptions,

        /// Apply the standard Tn5 offset.
        #[arg(long, action = clap::ArgAction::SetTrue)]
        tn5_shift: bool,
    },

    /// Infer scaling factor and library size from a normalised BigWig file.
    #[command(name = "bigwig-infer-scale", visible_alias = "infer-scale")]
    InferScale {
        /// Input BigWig file (CPM or RPKM normalised).
        #[arg(short, long, value_name = "BW")]
        bigwig: PathBuf,

        /// Output format.
        #[arg(long, value_enum, default_value = "table")]
        format: InferScaleFormat,

        /// Write output to FILE instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        /// Stop scanning a chromosome after N consecutive intervals with no new minimum.
        /// Set to 0 to disable early exit (full scan). Default 50000.
        #[arg(long, value_name = "N", default_value = "50000")]
        min_stable_streak: u32,

        /// Maximum interval start positions collected for bin-size (GCD) inference.
        /// Default 512 — GCD stabilises within the first few entries.
        #[arg(long, value_name = "N", default_value = "512")]
        max_starts: usize,
    },

    /// Compute between-sample scaling factors from BAM files.
    #[command(
        name = "bam-normalize",
        visible_alias = "norm-factors",
        after_help = "Example:\n  bamnado bam-normalize -b chip1.bam -b chip2.bam \\\n    --method csaw-background --bin-size 10000 --blacklist blacklist.bed"
    )]
    BamNormalize {
        /// Input BAM files. Repeat the flag for each sample.
        #[arg(short, long, value_name = "BAM")]
        bams: Vec<PathBuf>,

        /// CSV sample sheet with `sample` and `bam` columns. Alternative to --bams.
        #[arg(long, value_name = "CSV")]
        sample_sheet: Option<PathBuf>,

        /// Normalization method.
        #[arg(
            short,
            long,
            value_enum,
            default_value = "csaw-background",
            value_name = "METHOD"
        )]
        method: crate::normalization_factors::NormFactorMethod,

        /// Bin size for background windows.
        #[arg(short = 's', long, default_value = "10000", value_name = "BP")]
        bin_size: u64,

        /// Drop this percentage of highest-count bins before estimating (likely enriched).
        #[arg(long, default_value = "5.0", value_name = "PCT")]
        exclude_top_percent: f64,

        /// Sample name to use as the TMM reference. Defaults to the library size closest to the geometric mean.
        #[arg(long, value_name = "NAME")]
        reference_sample: Option<String>,

        /// TMM trim fraction for M-values.
        #[arg(long, default_value = "0.3", value_name = "FRACTION")]
        logratio_trim: f64,

        /// TMM trim fraction for A-values.
        #[arg(long, default_value = "0.05", value_name = "FRACTION")]
        sum_trim: f64,

        /// Reference-name prefix identifying spike-in sequences (--method spike-in).
        #[arg(short, long, value_name = "PREFIX")]
        exogenous_prefix: Option<String>,

        /// Count fragments instead of individual read alignments.
        #[arg(
            long = "fragment-counts",
            visible_alias = "use-fragment",
            action = clap::ArgAction::SetTrue
        )]
        use_fragment: bool,

        /// Skip scaffold or unplaced chromosomes.
        #[arg(
            long = "ignore-scaffolds",
            visible_alias = "ignore-scaffold",
            action = clap::ArgAction::SetTrue
        )]
        ignore_scaffold: bool,

        /// Write the bin count matrix to this TSV for QC (e.g. M-A plots).
        #[arg(long, value_name = "TSV")]
        counts_out: Option<PathBuf>,

        /// Output format.
        #[arg(long, value_enum, default_value = "table")]
        format: NormalizeFormat,

        /// Write output to FILE instead of stdout.
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,

        #[command(flatten)]
        filter_options: FilterOptions,
    },
}

#[derive(Debug, Clone, clap::ValueEnum)]
enum InferScaleFormat {
    Table,
    Tsv,
    Json,
}

#[derive(Debug, Clone, clap::ValueEnum)]
enum NormalizeFormat {
    Table,
    Tsv,
    Json,
}

/// Read a CSV sample sheet with `sample` and `bam` columns.
fn read_sample_sheet(path: &Path) -> Result<(Vec<PathBuf>, Vec<String>)> {
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.to_path_buf()))?
        .finish()
        .context("Failed to read sample sheet")?;

    let sample_col = df
        .column("sample")
        .context("Missing 'sample' column in sample sheet")?
        .str()
        .context("'sample' column is not string type")?;
    let bam_col = df
        .column("bam")
        .context("Missing 'bam' column in sample sheet")?
        .str()
        .context("'bam' column is not string type")?;

    let mut samples = Vec::new();
    let mut bams = Vec::new();
    for (sample, bam) in sample_col.into_iter().zip(bam_col) {
        let sample = sample.ok_or_else(|| anyhow::anyhow!("Null value in 'sample' column"))?;
        let bam = bam.ok_or_else(|| anyhow::anyhow!("Null value in 'bam' column"))?;
        samples.push(sample.to_string());
        bams.push(PathBuf::from(bam));
    }

    Ok((bams, samples))
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

fn log_active_filters(filter_options: &FilterOptions) {
    let mut rows = vec![
        ("strand", filter_options.strand.to_string()),
        ("min MAPQ", filter_options.min_mapq.to_string()),
        (
            "read length",
            format!(
                "{}..{} bp",
                filter_options.min_length, filter_options.max_length
            ),
        ),
        (
            "proper pairs",
            if filter_options.proper_pair {
                "required".to_string()
            } else {
                "off".to_string()
            },
        ),
    ];

    if filter_options.min_fragment_length.is_some() || filter_options.max_fragment_length.is_some()
    {
        rows.push((
            "fragment length",
            format!(
                "{}..{} bp",
                filter_options
                    .min_fragment_length
                    .map_or(String::from("0"), |v| v.to_string()),
                filter_options
                    .max_fragment_length
                    .map_or(String::from("inf"), |v| v.to_string())
            ),
        ));
    }
    if let Some(blacklist) = &filter_options.blacklisted_locations {
        rows.push(("blacklist", blacklist.display().to_string()));
    }
    if let Some(whitelist) = &filter_options.whitelisted_barcodes {
        rows.push(("barcode allowlist", whitelist.display().to_string()));
    }
    if let Some(rg) = &filter_options.read_group {
        rows.push(("read group", rg.clone()));
    }
    if let Some(tag) = &filter_options.filter_tag {
        let value = filter_options.filter_tag_value.as_deref().unwrap_or("*");
        rows.push(("tag", format!("{tag}={value}")));
    }
    if filter_options.ignore_duplicates {
        rows.push(("duplicates", "excluded".to_string()));
    }
    if filter_options.ignore_secondary {
        rows.push(("secondary alignments", "excluded".to_string()));
    }
    if filter_options.ignore_supplementary {
        rows.push(("supplementary alignments", "excluded".to_string()));
    }

    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL)
        .apply_modifier(UTF8_ROUND_CORNERS)
        .set_content_arrangement(ContentArrangement::Dynamic)
        .set_header(vec!["filter", "value"]);

    for (label, value) in rows {
        table.add_row(vec![
            Cell::new(label),
            Cell::new(value).set_alignment(CellAlignment::Right),
        ]);
    }

    info!("Read Filter Settings\n{table}");
}

fn create_filter_from_options(
    filter_options: &FilterOptions,
    bam_stats: Option<&crate::BamStats>,
) -> Result<crate::read_filter::BamReadFilter> {
    // Process whitelisted barcodes
    let whitelisted_barcodes = match &filter_options.whitelisted_barcodes {
        Some(whitelist) => {
            let barcodes =
                crate::CellBarcodes::from_csv(whitelist).context("Failed to read barcodes")?;
            Some(barcodes.barcodes())
        }
        None => None,
    };

    // Process blacklisted locations
    let blacklisted_locations = match (&filter_options.blacklisted_locations, bam_stats) {
        (Some(blacklist), Some(stats)) => {
            let lapper = crate::bam_utils::bed_to_lapper(blacklist.clone())
                .context("Failed to read blacklisted locations")?;

            let lapper = crate::bam_utils::convert_lapper_chrom_names_to_ids(lapper, stats)
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

    // Fragment length filtering only makes sense for paired-end data.
    if (filter_options.min_fragment_length.is_some()
        || filter_options.max_fragment_length.is_some())
        && let Some(stats) = bam_stats
        && !stats.is_paired_end()?
    {
        return Err(anyhow::anyhow!(
            "Fragment length filtering (--min-fragment-length / --max-fragment-length) \
                     requires paired-end reads, but the BAM file does not appear to contain \
                     paired-end data."
        ));
    }

    Ok(crate::read_filter::BamReadFilter::new(
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
        filter_options.min_fragment_length,
        filter_options.max_fragment_length,
    )
    .with_exclude_duplicates(filter_options.ignore_duplicates)
    .with_exclude_secondary(filter_options.ignore_secondary)
    .with_exclude_supplementary(filter_options.ignore_supplementary))
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

/// Parse `args` and run the CLI, returning the process exit code.
pub fn run_cli<I, T>(args: I) -> i32
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = match Cli::try_parse_from(args) {
        Ok(cli) => cli,
        Err(error) => {
            let _ = error.print();
            return error.exit_code();
        }
    };

    match dispatch(cli) {
        Ok(()) => 0,
        Err(error) => {
            eprintln!("Error: {error:?}");
            1
        }
    }
}

fn init_logging(verbose: u8) {
    let log_level = match verbose {
        0 => log::LevelFilter::Error,
        1 => log::LevelFilter::Warn,
        2 => log::LevelFilter::Info,
        3 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };

    let _ = colog::default_builder().filter_level(log_level).try_init();
    log::set_max_level(log_level);
}

fn dispatch(cli: Cli) -> Result<()> {
    init_logging(cli.verbose);

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
                crate::BamStats::new(bam.clone()).context("Failed to read BAM stats")?;

            // Log active filters
            log_active_filters(filter_options);

            // Create filter
            let filter = create_filter_from_options(filter_options, Some(&bam_stats))?;

            // Create pileup
            let coverage = crate::coverage_analysis::BamPileup::new(
                bam.clone(),
                coverage_options.bin_size.unwrap_or(50),
                coverage_options
                    .norm_method
                    .clone()
                    .unwrap_or(crate::signal_normalization::NormalizationMethod::Raw),
                coverage_options.scale_factor.unwrap_or(1.0),
                coverage_options.use_fragment,
                filter,
                true,
                coverage_options.ignore_scaffold,
                coverage_options.shift,
                coverage_options.truncate,
                Some(coverage_options.threads),
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

            // Log active filters
            log_active_filters(filter_options);

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
                    let barcodes = crate::bam_utils::CellBarcodesMulti::from_csv(whitelist)
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
                    crate::BamStats::new(bam.clone()).context("Failed to read BAM stats")?;

                // Process blacklisted locations
                let blacklisted_locations = match &filter_options.blacklisted_locations {
                    Some(blacklist) => {
                        let lapper = crate::bam_utils::bed_to_lapper(blacklist.clone())
                            .context("Failed to read blacklisted locations")?;
                        let lapper =
                            crate::bam_utils::convert_lapper_chrom_names_to_ids(lapper, &bam_stats)
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
                let filter = crate::read_filter::BamReadFilter::new(
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
                    filter_options.min_fragment_length,
                    filter_options.max_fragment_length,
                )
                .with_exclude_duplicates(filter_options.ignore_duplicates)
                .with_exclude_secondary(filter_options.ignore_secondary)
                .with_exclude_supplementary(filter_options.ignore_supplementary);

                // Create pileup
                pileups.push(crate::coverage_analysis::BamPileup::new(
                    bam.clone(),
                    coverage_options.bin_size.unwrap_or(50),
                    coverage_options
                        .norm_method
                        .clone()
                        .unwrap_or(crate::signal_normalization::NormalizationMethod::Raw),
                    coverage_options.scale_factor.unwrap_or(1.0),
                    coverage_options.use_fragment,
                    filter,
                    false,
                    coverage_options.ignore_scaffold,
                    coverage_options.shift,
                    coverage_options.truncate,
                    Some(coverage_options.threads),
                ));
            }

            // Create multi-BAM pileup
            let coverage = crate::coverage_analysis::MultiBamPileup::new(pileups);

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
            filter_options,
        } => {
            // Validate input BAM file
            validate_bam_file(input)?;

            // Log active filters
            log_active_filters(filter_options);

            // Create filter
            let bam_stats =
                crate::BamStats::new(input.clone()).context("Failed to read BAM stats")?;
            let filter = create_filter_from_options(filter_options, Some(&bam_stats))?;

            // Create and run BAM splitter
            let mut split = crate::spike_in_analysis::BamSplitter::new(
                input.clone(),
                output.clone(),
                exogenous_prefix.clone(),
                filter,
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

            // Log active filters
            log_active_filters(filter_options);

            // Create filter
            let bam_stats =
                crate::BamStats::new(input.clone()).context("Failed to read BAM stats")?;
            let filter = create_filter_from_options(filter_options, Some(&bam_stats))?;

            // Create BAM splitter
            let split =
                crate::bam_splitter::BamFilterer::new(input.clone(), output.clone(), filter);

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

            // Log active filters
            log_active_filters(filter_options);

            // Create filter
            let bam_stats =
                crate::BamStats::new(input.clone()).context("Failed to read BAM stats")?;
            let filter = create_filter_from_options(filter_options, Some(&bam_stats))?;

            // Create BAM modifier
            let modifier = crate::bam_modifier::BamModifier::new(
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
            pseudocount,
            scale_factor_bw1,
            scale_factor_bw2,
            threads,
        } => {
            crate::bigwig_compare::compare_bigwigs_with_threads(
                bw1,
                bw2,
                output,
                crate::bigwig_compare::CompareOptions {
                    comparison: comparison.clone(),
                    bin_size: *bin_size,
                    pseudocount: *pseudocount,
                    scale_factor1: *scale_factor_bw1,
                    scale_factor2: *scale_factor_bw2,
                },
                *threads,
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
            pseudocount,
            scale_factors,
            threads,
        } => {
            if bigwigs.is_empty() {
                return Err(anyhow::anyhow!("At least one BigWig file must be provided"));
            }

            let sfs = if scale_factors.is_empty() {
                None
            } else {
                Some(scale_factors.clone())
            };

            crate::bigwig_compare::aggregate_bigwigs_with_threads(
                bigwigs,
                output,
                method.clone(),
                *bin_size,
                *pseudocount,
                sfs,
                *threads,
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

            let collapsed = crate::bedgraph_utils::collapse_adjacent_bins(df)?;

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
        Commands::InferScale {
            bigwig,
            format,
            output,
            min_stable_streak,
            max_starts,
        } => {
            let config = crate::bigwig_infer_scale::InferScaleConfig {
                min_stable_streak: if *min_stable_streak == 0 {
                    u32::MAX
                } else {
                    *min_stable_streak
                },
                max_starts: *max_starts,
            };
            let result = crate::bigwig_infer_scale::infer_scale_factor(bigwig, &config)
                .context("Failed to infer scale factor from BigWig")?;

            // Build output string then write to file or stdout.
            let mut buf = String::new();
            match format {
                InferScaleFormat::Table => {
                    let mut table = Table::new();
                    table
                        .load_preset(UTF8_FULL)
                        .apply_modifier(UTF8_ROUND_CORNERS)
                        .set_content_arrangement(ContentArrangement::Dynamic)
                        .set_header(vec!["field", "value"]);

                    let rows: Vec<(&str, String)> = vec![
                        ("normalisation", result.norm_method.to_string()),
                        ("scale_factor", format!("{:.6e}", result.scale_factor)),
                        ("to_raw", format!("value × {:.6e}", result.scale_factor)),
                        ("library_size", format!("{:.0}", result.library_size)),
                        ("bin_size", format!("{} bp", result.bin_size)),
                        ("min_val", format!("{:.6e}", result.min_val)),
                        (
                            "second_min_val",
                            if result.second_min_val < f64::MAX {
                                format!("{:.6e}", result.second_min_val)
                            } else {
                                "n/a".to_string()
                            },
                        ),
                        (
                            "ratio",
                            if result.ratio.is_nan() {
                                "n/a".to_string()
                            } else {
                                format!("{:.4}", result.ratio)
                            },
                        ),
                        (
                            "pseudocount",
                            result
                                .pseudocount
                                .map(|p| format!("{:.4}", p))
                                .unwrap_or_else(|| "none".to_string()),
                        ),
                        (
                            "confident",
                            if result.confident { "yes" } else { "no" }.to_string(),
                        ),
                        ("chroms_scanned", result.chroms_scanned.to_string()),
                    ];

                    for (label, value) in &rows {
                        table.add_row(vec![
                            Cell::new(label),
                            Cell::new(value).set_alignment(CellAlignment::Right),
                        ]);
                    }

                    use std::fmt::Write as _;
                    writeln!(buf, "{table}").unwrap();

                    // Warnings always go to stderr regardless of --output.
                    for w in result.warnings.messages() {
                        eprintln!("Warning: {w}");
                    }
                }
                InferScaleFormat::Tsv => {
                    use std::fmt::Write as _;
                    writeln!(buf, "field\tvalue").unwrap();
                    writeln!(buf, "normalisation\t{}", result.norm_method).unwrap();
                    writeln!(buf, "scale_factor\t{:.6e}", result.scale_factor).unwrap();
                    writeln!(buf, "to_raw\tvalue × {:.6e}", result.scale_factor).unwrap();
                    writeln!(buf, "library_size\t{:.0}", result.library_size).unwrap();
                    writeln!(buf, "bin_size\t{}", result.bin_size).unwrap();
                    writeln!(buf, "min_val\t{:.6e}", result.min_val).unwrap();
                    writeln!(
                        buf,
                        "second_min_val\t{}",
                        if result.second_min_val < f64::MAX {
                            format!("{:.6e}", result.second_min_val)
                        } else {
                            "n/a".to_string()
                        }
                    )
                    .unwrap();
                    writeln!(
                        buf,
                        "ratio\t{}",
                        if result.ratio.is_nan() {
                            "n/a".to_string()
                        } else {
                            format!("{:.4}", result.ratio)
                        }
                    )
                    .unwrap();
                    writeln!(
                        buf,
                        "pseudocount\t{}",
                        result
                            .pseudocount
                            .map(|p| format!("{:.4}", p))
                            .unwrap_or_else(|| "none".to_string())
                    )
                    .unwrap();
                    writeln!(
                        buf,
                        "confident\t{}",
                        if result.confident { "yes" } else { "no" }
                    )
                    .unwrap();
                    writeln!(buf, "chroms_scanned\t{}", result.chroms_scanned).unwrap();
                    for w in result.warnings.messages() {
                        writeln!(buf, "warning\t{w}").unwrap();
                    }
                }
                InferScaleFormat::Json => {
                    buf = serde_json::to_string_pretty(&result)?;
                    buf.push('\n');
                }
            }

            match output {
                Some(path) => std::fs::write(path, &buf)
                    .with_context(|| format!("Cannot write to {}", path.display()))?,
                None => print!("{buf}"),
            }
        }

        Commands::BamNormalize {
            bams,
            sample_sheet,
            method,
            bin_size,
            exclude_top_percent,
            reference_sample,
            logratio_trim,
            sum_trim,
            exogenous_prefix,
            use_fragment,
            ignore_scaffold,
            counts_out,
            format,
            output,
            filter_options,
        } => {
            use crate::normalization_factors::NormFactorMethod;

            let (bam_paths, sample_names) = match (sample_sheet, bams.is_empty()) {
                (Some(sheet), true) => read_sample_sheet(sheet)?,
                (None, false) => {
                    let mut names: Vec<String> = bams
                        .iter()
                        .map(|b| {
                            b.file_stem()
                                .and_then(|s| s.to_str())
                                .map(String::from)
                                .unwrap_or_else(|| b.display().to_string())
                        })
                        .collect();

                    // Disambiguate duplicate sample names (e.g. the same BAM passed twice)
                    // so per-sample columns/labels stay unique downstream.
                    let mut occurrences: std::collections::HashMap<String, usize> =
                        std::collections::HashMap::new();
                    for name in &names {
                        *occurrences.entry(name.clone()).or_insert(0) += 1;
                    }
                    if occurrences.values().any(|&count| count > 1) {
                        let mut seen_so_far: std::collections::HashMap<String, usize> =
                            std::collections::HashMap::new();
                        for name in names.iter_mut() {
                            let total = occurrences[name.as_str()];
                            let seen = seen_so_far.entry(name.clone()).or_insert(0);
                            *seen += 1;
                            if total > 1 {
                                *name = format!("{name}_{seen}");
                            }
                        }
                    }

                    (bams.clone(), names)
                }
                (Some(_), false) => {
                    return Err(anyhow::anyhow!(
                        "Provide either --bams or --sample-sheet, not both"
                    ));
                }
                (None, true) => {
                    return Err(anyhow::anyhow!("Provide either --bams or --sample-sheet"));
                }
            };

            let min_samples = match method {
                NormFactorMethod::Tmm
                | NormFactorMethod::CsawBackground
                | NormFactorMethod::MedianOfRatios => 2,
                NormFactorMethod::Cpm | NormFactorMethod::SpikeIn => 1,
            };
            anyhow::ensure!(
                bam_paths.len() >= min_samples,
                "{method:?} requires at least {min_samples} sample(s), got {}",
                bam_paths.len()
            );

            for bam in &bam_paths {
                validate_bam_file(bam)?;
            }

            log_active_filters(filter_options);

            let bam_stats: Vec<crate::BamStats> = bam_paths
                .iter()
                .map(|p| {
                    crate::BamStats::new(p.clone())
                        .with_context(|| format!("Failed to read BAM stats for {}", p.display()))
                })
                .collect::<Result<Vec<_>>>()?;

            let filters: Vec<crate::read_filter::BamReadFilter> = bam_stats
                .iter()
                .map(|stats| create_filter_from_options(filter_options, Some(stats)))
                .collect::<Result<Vec<_>>>()?;

            let result = if matches!(method, NormFactorMethod::SpikeIn) {
                let prefix = exogenous_prefix.clone().ok_or_else(|| {
                    anyhow::anyhow!("--exogenous-prefix is required for --method spike-in")
                })?;
                crate::normalization_factors::spike_in_scale_factors(
                    &bam_stats,
                    &sample_names,
                    &prefix,
                )
                .context("Failed to compute spike-in scale factors")?
            } else {
                let counts = crate::bin_counts::count_bins(
                    &bam_paths,
                    &sample_names,
                    *bin_size,
                    &filters,
                    *use_fragment,
                    *ignore_scaffold,
                )
                .context("Failed to count reads per bin")?;

                if let Some(counts_out_path) = counts_out {
                    let mut df = counts.to_polars()?;
                    let mut file = std::fs::File::create(counts_out_path)
                        .context("Failed to create --counts-out file")?;
                    CsvWriter::new(&mut file)
                        .include_header(true)
                        .with_separator(b'\t')
                        .finish(&mut df)?;
                }

                let params = crate::normalization_factors::NormFactorParams {
                    reference_sample: reference_sample.clone(),
                    logratio_trim: *logratio_trim,
                    sum_trim: *sum_trim,
                    exclude_top_percent: *exclude_top_percent,
                    ..Default::default()
                };

                crate::normalization_factors::compute_scale_factors(&counts, method, &params)
                    .context("Failed to compute scale factors")?
            };

            let mut buf = String::new();
            match format {
                NormalizeFormat::Table => {
                    let mut table = Table::new();
                    table
                        .load_preset(UTF8_FULL)
                        .apply_modifier(UTF8_ROUND_CORNERS)
                        .set_content_arrangement(ContentArrangement::Dynamic)
                        .set_header(vec![
                            "sample",
                            "library_size",
                            "norm_factor",
                            "scale_factor",
                            "reference",
                        ]);

                    for i in 0..result.sample_names.len() {
                        let is_ref = result.reference_sample.as_deref()
                            == Some(result.sample_names[i].as_str());
                        table.add_row(vec![
                            Cell::new(&result.sample_names[i]),
                            Cell::new(result.library_sizes[i].to_string())
                                .set_alignment(CellAlignment::Right),
                            Cell::new(format!("{:.6}", result.norm_factors[i]))
                                .set_alignment(CellAlignment::Right),
                            Cell::new(format!("{:.6}", result.scale_factors[i]))
                                .set_alignment(CellAlignment::Right),
                            Cell::new(if is_ref { "yes" } else { "-" })
                                .set_alignment(CellAlignment::Right),
                        ]);
                    }

                    use std::fmt::Write as _;
                    writeln!(buf, "{table}").unwrap();
                    writeln!(buf, "method: {}", result.method).unwrap();
                    writeln!(
                        buf,
                        "bins: {} used / {} total",
                        result.n_bins_used, result.n_bins_total
                    )
                    .unwrap();
                }
                NormalizeFormat::Tsv => {
                    use std::fmt::Write as _;
                    writeln!(
                        buf,
                        "sample\tlibrary_size\tnorm_factor\tscale_factor\treference"
                    )
                    .unwrap();
                    for i in 0..result.sample_names.len() {
                        let is_ref = result.reference_sample.as_deref()
                            == Some(result.sample_names[i].as_str());
                        writeln!(
                            buf,
                            "{}\t{}\t{:.6}\t{:.6}\t{}",
                            result.sample_names[i],
                            result.library_sizes[i],
                            result.norm_factors[i],
                            result.scale_factors[i],
                            if is_ref { "yes" } else { "-" }
                        )
                        .unwrap();
                    }
                }
                NormalizeFormat::Json => {
                    buf = serde_json::to_string_pretty(&result)?;
                    buf.push('\n');
                }
            }

            match output {
                Some(path) => std::fs::write(path, &buf)
                    .with_context(|| format!("Cannot write to {}", path.display()))?,
                None => print!("{buf}"),
            }
        }
    }

    Ok(())
}
