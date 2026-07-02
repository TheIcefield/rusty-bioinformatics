use clap::Parser;

use fna::FnaFile;
use nucleotides::sequence::Sequence;
use plot_collection::PlotCpgIslandsOptions;
use plotters::prelude::*;

use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to FNA file
    #[arg(short, long)]
    input: PathBuf,

    /// Dir of output chars
    #[arg(short, long)]
    plot: Option<PathBuf>,

    /// Minimum length of CpG islands in nucleotides
    #[arg(long, default_value_t = 200)]
    cgp_min_len: usize,

    /// Minimum length of CpG islands in nucleotides
    #[arg(long, default_value_t = 10)]
    cgp_step: usize,

    /// Minimum GC content of CpG islands
    #[arg(long, default_value_t = 50.0)]
    cpg_min_gc: f64,

    /// Minimum Observed/Expected CG content of CpG islands
    #[arg(long, default_value_t = 0.6)]
    cpg_min_oe: f64,

    /// Minimum length of ORF in nucleotides
    #[arg(long, default_value_t = 300)]
    orf_min_len: usize,

    /// If ORF#1 is inside of ORF#2 keep only ORF#2
    #[arg(long, default_value_t = false)]
    best_orf_only: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    let fna = FnaFile::open(&cli.input)?;

    if let Some(charts_dir) = &cli.plot {
        let _ = std::fs::create_dir(charts_dir);
    }

    fna.records
        .iter()
        .enumerate()
        .map(|(record_id, record)| {
            let islands = record.content.find_cpg_islands(
                cli.cgp_min_len,
                cli.cgp_step,
                cli.cpg_min_gc,
                cli.cpg_min_oe,
            );

            let orfs = record.content.find_orfs(cli.orf_min_len);
            let orfs = if cli.best_orf_only { Sequence::best_orfs(orfs) } else { orfs };

            let intersections = record.content.find_subsequence_intersections(&islands, &orfs);

            let metrics = record.content.get_window_gc_oe_metrics(std::cmp::min(cli.cgp_min_len, cli.orf_min_len), cli.cgp_step);

            (record_id, record, islands, orfs, intersections, metrics)
        })
        .try_for_each(
            |(record_id, record, islands, orfs, intersections, metrics)| -> Result<(), Box<dyn std::error::Error>> {
                println!("#{record_id}: {}.", record.header);
                println!("GC: {}%", record.content.get_gc_content());

                println!("\nCpG-islands:");
                islands.iter().enumerate().for_each(|(island_id, island)| {
                    println!(
                        "    #{island_id}, start: {}, end: {}, len: {}, GC: {}%",
                        island.start,
                        island.end,
                        island.seq.length(),
                        island.seq.get_gc_content()
                    )
                });

                println!("\nORFs:");
                orfs.iter().enumerate().for_each(|(orf_id, orf)| {
                    println!(
                        "    #{orf_id}, start: {}, end: {}, len: {}, GC: {}%",
                        orf.start,
                        orf.end,
                        orf.seq.length(),
                        orf.seq.get_gc_content()
                    )
                });

                println!("\nIntersections:");
                intersections.iter().enumerate().for_each(|(intersection_id, (intersection, island_id, orf_id))| {
                    println!(
                        "    #{intersection_id}, CpG-Island: #{island_id}, ORF: #{orf_id}, start: {}, end: {}, len: {}, GC: {}%",
                        intersection.start,
                        intersection.end,
                        intersection.seq.length(),
                        intersection.seq.get_gc_content()
                    )
                });
                println!();

                if let Some(charts_dir) = &cli.plot {
                    let chart_path = charts_dir.join(format!("{record_id}.png"));
                    let root =
                        BitMapBackend::new(&chart_path, (1500, 1000)).into_drawing_area();

                    root.fill(&WHITE)?;

                    let opts = PlotCpgIslandsOptions {
                        show_gc: true,
                        show_oe: true,
                        show_min_gc: true,
                        show_min_oe: true,
                        show_cpg_islands: true,
                        show_orfs: true,
                        min_gc: cli.cpg_min_gc,
                        min_oe: cli.cpg_min_oe,
                    };

                    let title =  format!("{} #{record_id}", record.header);

                    plot_collection::plot_cpg_and_orf(
                        &root,
                        &title,
                        &metrics,
                        &Some(&islands),
                        &Some(&orfs),
                        &opts
                    )?;

                    root.present()?;
                }

                Ok(())
            },
        )?;

    Ok(())
}
