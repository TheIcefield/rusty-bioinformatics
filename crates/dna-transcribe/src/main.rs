use std::{io::Write, path::PathBuf};

use clap::Parser;
use fna::{FnaFile, FnaRecord};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to FNA file
    #[arg(short, long)]
    input: PathBuf,

    /// Path to transcribed FNA file
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let mut output = std::fs::File::create(
        cli.output
            .unwrap_or(cli.input.with_extension("transcribed.fna")),
    )?;

    let fna = FnaFile::open(&cli.input)?;

    for record in fna.records {
        let transcribed_record = FnaRecord {
            header: record.header,
            content: record.content.transcribe()?,
        };

        writeln!(&mut output, "{}", transcribed_record)?;
    }

    Ok(())
}
