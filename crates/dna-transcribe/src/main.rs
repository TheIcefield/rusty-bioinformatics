use std::{fs::File, io::Write, path::PathBuf};

use clap::Parser;
use fna::{FnaFile, FnaRecord};
use nucleotides::sequence::Sequence;

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

fn write_seq(file: &mut File, seq: &Sequence) -> Result<(), Box<dyn std::error::Error>> {
    const CHUNK_SIZE: usize = 60;

    let mut it = seq.0.iter();
    let mut cur_pos = 0usize;

    while cur_pos < seq.length() {
        let start = cur_pos;
        let mut count = 0;

        while count < CHUNK_SIZE {
            let Some(_) = it.next() else {
                break;
            };

            count += 1;
            cur_pos += 1;
        }

        let sub_seq = &seq.0[start..cur_pos];
        let sub_str = Sequence::to_string(sub_seq);

        writeln!(file, "{sub_str}")?;
    }

    Ok(())
}

fn write_record(file: &mut File, record: &FnaRecord) -> Result<(), Box<dyn std::error::Error>> {
    writeln!(file, ">{}", record.header)?;
    write_seq(file, &record.content)
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

        write_record(&mut output, &transcribed_record)?;
    }

    Ok(())
}
