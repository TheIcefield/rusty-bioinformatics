use std::path::PathBuf;

use clap::Args;
use fna::FnaFile;

#[derive(Args, Debug)]
pub struct CompareSubcommand {
    #[arg(long)]
    second: PathBuf,

    #[arg(long)]
    hamming_distance: bool,

    #[arg(long)]
    transitions: bool,

    #[arg(long)]
    transversions: bool,
}

pub fn process_compare_cmd(
    fna: &FnaFile,
    args: &CompareSubcommand,
) -> Result<(), Box<dyn std::error::Error>> {
    let second_fna = FnaFile::open(&args.second)?;

    if fna.records.len() != second_fna.records.len() {
        return Err(format!(
            "Records count is different: {} and {}",
            fna.records.len(),
            second_fna.records.len()
        )
        .into());
    }

    fna.records
        .iter()
        .zip(second_fna.records.iter())
        .enumerate()
        .try_for_each(
            |(record_id, (first, second))| -> Result<(), Box<dyn std::error::Error>> {
                println!(
                    "Compare record #{record_id}\n    {}\nVS\n    {}\n",
                    first.header, second.header
                );

                if args.hamming_distance {
                    let distance = first.content.get_hamming_distance(&second.content)?;

                    println!("Hamming distance: {distance}");
                }

                let transitions_count = if args.transitions {
                    let transitions = first.content.get_transitions(&second.content)?;
                    println!("Transitions: {transitions}");
                    transitions
                } else {
                    0
                };

                let transversions_count = if args.transversions {
                    let transversions = first.content.get_transversions(&second.content)?;
                    println!("Transversions: {transversions}");
                    transversions
                } else {
                    0
                };

                if args.transitions && args.transversions && transversions_count != 0 {
                    let ratio = transitions_count as f64 / transversions_count as f64;
                    println!("Transitions/Transversions ratio: {ratio}");
                }

                Ok(())
            },
        )
}
