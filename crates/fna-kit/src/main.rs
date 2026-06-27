use fna::{FnaFile, FnaRecord};

use clap::{Parser, Subcommand};

use std::{
    io::Write,
    path::{Path, PathBuf},
};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: ActionSubcommand,
}

#[derive(Subcommand, Debug)]
enum ActionSubcommand {
    /// Split records of FNA file into separate FNA files.
    SplitRecords {
        /// Path to FNA file.
        #[arg(short, long)]
        input: PathBuf,
    },
    /// Merge record of specified FNA file into one record written in FNA file.
    MergeRecords {
        /// Path to FNA file.
        #[arg(short, long)]
        input: PathBuf,

        /// Path to output FNA file.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Merge specified files into one FNA file.
    MergeFiles {
        /// Path to FNA files.
        #[arg(short, long)]
        input: Vec<PathBuf>,

        /// Path to output FNA file.
        #[arg(short, long)]
        output: PathBuf,
    },
}

fn split_records(input: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let fna = FnaFile::open(input)?;

    for (record_id, record) in fna.records.into_iter().enumerate() {
        let file_name = input.with_extension(format!("{record_id}.fna"));
        let mut out = std::fs::File::create(&file_name)?;

        writeln!(out, "{record}")?;
    }

    Ok(())
}

fn merge_records(input: &Path, output: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let fna = FnaFile::open(input)?;

    let mut full_record = FnaRecord::default();

    for mut record in fna.records.into_iter() {
        full_record.header.push_str(&format!("{}. ", record.header));
        full_record.content.0.append(&mut record.content.0);
    }

    let mut out = std::fs::File::create(output)?;
    writeln!(out, "{full_record}")?;

    Ok(())
}

fn merge_files(input: &[PathBuf], output: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let mut records = Vec::<FnaRecord>::new();

    for input in input.iter() {
        let fna = FnaFile::open(input)?;

        for record in fna.records {
            records.push(record);
        }
    }

    let mut out = std::fs::File::create(output)?;
    for record in records {
        writeln!(out, "{record}")?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        ActionSubcommand::SplitRecords { input } => split_records(input)?,
        ActionSubcommand::MergeRecords { input, output } => merge_records(
            input,
            &output.clone().unwrap_or(input.with_extension("merged.fna")),
        )?,
        ActionSubcommand::MergeFiles { input, output } => merge_files(input, output)?,
    }

    Ok(())
}
