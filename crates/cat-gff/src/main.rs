use gff::GffFile;

use clap::Parser;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to GFF file
    #[arg(short, long)]
    file: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let gff = GffFile::open(&cli.file)?;

    gff.features
        .iter()
        .enumerate()
        .for_each(|(feature_id, feature)| {
            println!("Feature #{feature_id}");
            println!("  SeqId: {}", feature.seq_id);
            println!("  Source: {}", feature.source);
            println!("  Feature type: {}", feature.feature_type);
            println!("  Position: {}-{}", feature.start, feature.end);
            println!("  Score: {}", feature.score);
            println!("  Strand: {}", feature.strand);
            println!("  Phase: {}", feature.phase);

            feature
                .attributes
                .iter()
                .enumerate()
                .for_each(|(attribute_id, (key, value))| {
                    println!("  Attribute #{attribute_id}: {key}={value}")
                });
        });

    Ok(())
}
