use fna::FnaFile;

use clap::{Parser, Subcommand};
use gff::GffFile;

use std::path::PathBuf;

pub mod compare_cmd;
pub mod consts;
pub mod find_cmd;
pub mod info_cmd;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to FNA file
    #[arg(short, long)]
    fna: PathBuf,

    /// Path to GFF file
    #[arg(short, long)]
    gff: Option<PathBuf>,

    #[command(subcommand)]
    command: ContentSubcommand,
}

#[derive(Subcommand, Debug)]
enum ContentSubcommand {
    /// Show speciefied info about given FNA
    Info(info_cmd::InfoSubcommandArgs),

    /// Compare two given FNA files
    Compare(compare_cmd::CompareSubcommand),

    /// Find special regions in FNA file
    #[command(subcommand)]
    Find(find_cmd::FindSubcommand),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let fna = FnaFile::open(&cli.fna)?;

    let gff = if let Some(p) = &cli.gff {
        Some(GffFile::open(p)?)
    } else {
        None
    };

    match &cli.command {
        ContentSubcommand::Info(args) => info_cmd::process_info_cmd(&fna, args),
        ContentSubcommand::Find(args) => find_cmd::process_find_cmd(&fna, &gff, args),
        ContentSubcommand::Compare(args) => compare_cmd::process_compare_cmd(&fna, args),
    }
}
