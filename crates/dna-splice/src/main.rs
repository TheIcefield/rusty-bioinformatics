use std::{fs::File, io::Write, path::PathBuf};

use clap::Parser;
use fna::FnaFile;
use nucleotides::{
    codons::{Codon, CodonSequence},
    sequence::Sequence,
};

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

fn write_seq(file: &mut File, seq: &[Codon]) -> Result<(), Box<dyn std::error::Error>> {
    const CHUNK_SIZE: usize = 60;

    let mut it = seq.iter();
    let mut cur_pos = 0usize;

    while cur_pos < seq.len() {
        let start = cur_pos;
        let mut count = 0;

        while count < CHUNK_SIZE {
            let Some(codon) = it.next() else {
                break;
            };

            cur_pos += 1;

            if *codon == Codon::Stop {
                count += 4;
            } else {
                count += 1;
            }
        }

        let Some(sub_seq) = seq.get(start..cur_pos) else {
            continue;
        };

        let sub_str = CodonSequence::to_string(sub_seq);

        writeln!(file, "{sub_str}")?;
    }

    Ok(())
}

fn write_record(
    file: &mut File,
    header: &str,
    content: &[Codon],
) -> Result<(), Box<dyn std::error::Error>> {
    writeln!(file, ">{}", header)?;
    write_seq(file, content)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let mut output = std::fs::File::create(
        cli.output
            .unwrap_or(cli.input.with_extension("spliced.fna")),
    )?;

    let fna = FnaFile::open(&cli.input)?;
    assert!(fna.records.len() > 1);

    let initial_seq = &fna.records[0];

    let introns = fna
        .records
        .iter()
        .skip(1)
        .map(|record| record.content.clone())
        .collect::<Vec<Sequence>>();

    println!("Registered introns... {}", introns.len());

    let spliced = initial_seq.content.splice(&introns);
    let transcribed = spliced.transcribe()?;
    let translated = transcribed.translate()?;

    write_record(&mut output, &initial_seq.header, &translated.0.0)?;

    Ok(())
}
