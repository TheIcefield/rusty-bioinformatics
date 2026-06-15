use fna::FnaFile;
use nucleotides::nucleotide::Nucleotide;

use clap::{Parser, Subcommand};
use plotters::{backend::BitMapBackend, prelude::*, style::full_palette::WHITE};

use std::{
    io::{Read, Write},
    path::{Path, PathBuf},
};

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
    Transcribe {
        #[arg(short, long)]
        out: PathBuf,
    },

    #[command(subcommand)]
    Find(find_cmd::FindSubcommand),
}

fn read_fna(path: &Path) -> Result<FnaFile, String> {
    let mut file = std::fs::File::open(path).unwrap();

    let mut raw_data = String::new();
    file.read_to_string(&mut raw_data)
        .map_err(|err| err.to_string())?;

    FnaFile::read(&raw_data)
}

fn main() {
    let cli = Cli::parse();

    let fna = match read_fna(&cli.file) {
        Ok(fna) => fna,
        Err(err) => {
            eprintln!("{err}");
            return;
        }
    };

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

        Some(ContentSubcommand::Transcribe { out }) => {
            let mut file = std::fs::File::create(out).unwrap();

            for record in fna.records.iter() {
                writeln!(file, "> {}", record.header).unwrap();

                match record.content.transcribe() {
                    Ok(seq) => writeln!(file, "{seq}").unwrap(),
                    Err(err) => eprintln!("{err}"),
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

                    let metrics = record.content.cpg_islands_metrics(*min_len, *step);

                    (record_idx, record, islands, metrics)
                })
                .for_each(
                    |(record_idx, record, islands, (positions, gc_values, oe_values))| {
                        println!("#{record_idx}: {}.", record.header);

                        for (island_idx, island) in islands.iter().enumerate() {
                            println!(
                                "    Island #{island_idx}, start: {}, end: {}",
                                island.start, island.end
                            );
                        }

                        if let Some(charts_dir) = plot {
                            let chart_path = charts_dir.join(format!("{record_idx}.png"));
                            let root =
                                BitMapBackend::new(&chart_path, (1200, 800)).into_drawing_area();

                            if let Err(err) = root.fill(&WHITE) {
                                eprintln!("{err}");
                                return;
                            }

                            if let Err(err) = find_cmd::plot_cpg_islands(
                                &root,
                                *min_gc,
                                *min_oe,
                                (&positions, &gc_values, &oe_values),
                                &islands,
                            ) {
                                eprintln!("{err}");
                            }
                        }
                    },
                );
        }

        None => eprintln!("Subcommand not provided!"),
    }
}
