use gff::GffFile;

use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to GFF file
    #[arg(short, long)]
    file: PathBuf,

    /// Displays range of features in range first..last (by default: 0..N)
    #[arg(long)]
    first: Option<usize>,

    /// Displays range of features in range first..last (by default: 0..N)
    #[arg(long)]
    last: Option<usize>,

    /// Displays seq-id
    #[arg(long)]
    seq_id: bool,

    /// Displays source
    #[arg(long)]
    source: bool,

    /// Displays feature-type
    #[arg(long)]
    feature_type: bool,

    /// Displays position and length
    #[arg(long)]
    position: bool,

    /// Displays score
    #[arg(long)]
    score: bool,

    /// Displays strand
    #[arg(long)]
    strand: bool,

    /// Displays phase
    #[arg(long)]
    phase: bool,

    /// Displays attributes
    #[arg(long)]
    attributes: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let gff = GffFile::open(&cli.file)?;

    gff.features
        .iter()
        .enumerate()
        .filter(|(id, _)| {
            *id >= cli.first.unwrap_or(0) && *id <= cli.last.unwrap_or(gff.features.len())
        })
        .for_each(|(feature_id, feature)| {
            println!("Feature #{feature_id}:");

            if cli.seq_id {
                println!("  SeqId: {}", feature.seq_id);
            }

            if cli.source {
                println!("  Source: {}", feature.source);
            }

            if cli.feature_type {
                println!("  Feature type: {}", feature.feature_type);
            }

            if cli.position {
                println!(
                    "  Position: {}-{}. (Total len: {})",
                    feature.start,
                    feature.end,
                    feature.end - feature.start
                );
            }

            if cli.score {
                println!("  Score: {}", feature.score);
            }

            if cli.strand {
                println!("  Strand: {}", feature.strand);
            }

            if cli.phase {
                println!("  Phase: {}", feature.phase);
            }

            if cli.attributes {
                feature
                    .attributes
                    .iter()
                    .enumerate()
                    .for_each(|(attribute_id, (key, value))| {
                        println!("  Attribute #{attribute_id}: {key}={value}")
                    });
            }
        });

    Ok(())
}
