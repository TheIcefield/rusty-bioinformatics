use fna::FnaFile;
use nucleotides::nucleotide::Nucleotide;

use clap::{Parser, Subcommand};
use plotters::{backend::BitMapBackend, prelude::*, style::full_palette::WHITE};

use std::path::PathBuf;

pub mod find_cmd;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to FNA file
    #[arg(short, long)]
    file: PathBuf,

    #[command(subcommand)]
    command: Option<ContentSubcommand>,
}

#[derive(Subcommand, Debug)]
enum ContentSubcommand {
    Info {
        /// Displays loaded sequence kind (DNA or RNA)
        #[arg(long)]
        kind: bool,

        /// Displays whole number of nucleotides
        #[arg(long)]
        len: bool,

        /// Displays number of each nucleotide
        #[arg(long)]
        count: bool,

        /// Displays GC content
        #[arg(long)]
        gc: bool,
    },

    #[command(subcommand)]
    Find(find_cmd::FindSubcommand),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let fna = FnaFile::open(&cli.file)?;

    match &cli.command {
        Some(ContentSubcommand::Info {
            kind,
            len,
            count,
            gc,
        }) => {
            for (idx, record) in fna.records.iter().enumerate() {
                println!("#{idx}: {}.", record.header);

                if *kind {
                    println!("    Sequence kind: {}", record.content.get_kind());
                }

                if *len {
                    println!("    Length: {} nucleotides", record.content.length());
                }

                if *count {
                    println!("    A: {}", record.content.count(Nucleotide::Adenine));

                    if record.content.is_dna() {
                        println!("    T: {}", record.content.count(Nucleotide::Thymine));
                    } else {
                        println!("    T: {}", record.content.count(Nucleotide::Uracil));
                    }

                    println!("    G: {}", record.content.count(Nucleotide::Guanine));
                    println!("    C: {}", record.content.count(Nucleotide::Cytosine));
                }

                if *gc {
                    println!("    GC content: {}%", record.content.get_gc_content());
                }
            }
        }

        Some(ContentSubcommand::Find(find_cmd::FindSubcommand::CpgIslands {
            min_len,
            min_gc,
            min_oe,
            step,
            plot,
        })) => {
            if let Some(charts_dir) = plot {
                let _ = std::fs::create_dir(charts_dir);
            }

            fna.records
                .iter()
                .enumerate()
                .map(|(record_idx, record)| {
                    let islands = record
                        .content
                        .find_cpg_islands(*min_len, *min_gc, *min_oe, *step);

                    let metrics = record.content.get_window_gc_oe_metrics(*min_len, *step);

                    (record_idx, record, islands, metrics)
                })
                .try_for_each(
                    |(record_idx, record, islands, (positions, gc_values, oe_values))| -> Result<(), Box<dyn std::error::Error>> {
                        println!("#{record_idx}: {}.", record.header);

                        for (island_id, island) in islands.iter().enumerate() {
                            println!(
                                "    Island #{island_id}, start: {}, end: {}, len: {}, GC: {}%",
                                island.start,
                                island.end,
                                island.seq.length(),
                                island.seq.get_gc_content()
                            );
                        }

                        if let Some(charts_dir) = plot {
                            let chart_path = charts_dir.join(format!("{record_idx}.png"));
                            let root =
                                BitMapBackend::new(&chart_path, (1200, 1000)).into_drawing_area();

                            root.fill(&WHITE)?;

                            find_cmd::plot_cpg_islands(
                                &root,
                                &format!("{} #{record_idx}", record.header),
                                *min_gc,
                                *min_oe,
                                (&positions, &gc_values, &oe_values),
                                &islands,
                            )?;

                            root.present()?;
                        }

                        Ok(())
                    },
                )?;
        }

        Some(ContentSubcommand::Find(find_cmd::FindSubcommand::Orf { min_len })) => {
            fna.records
                .iter()
                .enumerate()
                .filter_map(|(record_idx, record)| {
                    let orfs = record.content.find_orfs(*min_len);

                    match orfs {
                        Ok(orfs) => Some((record_idx, record, orfs)),
                        Err(_) => None,
                    }
                })
                .try_for_each(
                    |(record_idx, record, orfs)| -> Result<(), Box<dyn std::error::Error>> {
                        println!("#{record_idx}: {}.", record.header);

                        for (orf_id, orf) in orfs.iter().enumerate() {
                            println!(
                                "    ORF #{orf_id}, start: {}, end: {}, len: {}, GC: {}%",
                                orf.start,
                                orf.end,
                                orf.seq.length(),
                                orf.seq.get_gc_content()
                            );
                        }

                        Ok(())
                    },
                )?;
        }

        None => eprintln!("Subcommand not provided!"),
    }

    Ok(())
}
