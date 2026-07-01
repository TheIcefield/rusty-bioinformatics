use fna::FnaFile;
use nucleotides::sequence::{CpgIsland, Sequence};

use clap::{Args, Subcommand};
use plotters::{coord::Shift, prelude::*};

use std::{path::PathBuf, str::FromStr};

use crate::consts;

#[derive(Subcommand, Debug)]
pub enum FindSubcommand {
    CpgIslands(FindCpgsArgs),
    Orfs(FindOrfsArgs),
    Motifs(FindMotifsArgs),
}

#[derive(Args, Debug)]
pub struct FindCpgsArgs {
    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_LEN)]
    min_len: usize,

    #[arg(long, default_value_t = consts::DEFAULT_CPG_WIN_STEP)]
    step: usize,

    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_GC)]
    min_gc: f64,

    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_OE)]
    min_oe: f64,

    #[arg(long)]
    plot: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub struct FindOrfsArgs {
    #[arg(long, default_value_t = 300)]
    min_len: usize,
}

#[derive(Args, Debug)]
pub struct FindMotifsArgs {
    #[arg(long)]
    sub_str: String,
}

fn plot_cpg_islands(
    root: &DrawingArea<BitMapBackend, Shift>,
    title: &str,
    min_gc: f64,
    min_oe: f64,
    metrics: (&[usize], &[f64], &[f64]),
    islands: &[CpgIsland],
) -> Result<(), Box<dyn std::error::Error>> {
    const FONT_NAME: &str = "sans-serif";
    const TITLE_FONT_SIZE: f64 = 30f64;
    const CAPTION_FONT_SIZE: f64 = 20f64;
    const OE_MIN: f64 = 0.0;
    const OE_MAX: f64 = 2.0;
    const GC_MIN: f64 = 0.0;
    const GC_MAX: f64 = 100.0;

    let (positions, gc_values, oe_values) = metrics;
    let x_min = *positions.first().unwrap_or(&0) as f64;
    let x_max = *positions.last().unwrap_or(&0) as f64;

    root.titled(title, (FONT_NAME, TITLE_FONT_SIZE))?;

    let (upper, lower) = root.split_vertically(500);

    let mut oe_chart = ChartBuilder::on(&upper)
        .caption("Observed/Expected CpG", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .margin_top(30)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, OE_MIN..OE_MAX)?;

    let mut gc_chart = ChartBuilder::on(&lower)
        .caption("GC content (%)", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, GC_MIN..GC_MAX)?;

    plot_collection::plot_oe_values(&mut oe_chart, positions, oe_values, min_oe)?;
    plot_collection::plot_gc_values(&mut gc_chart, positions, gc_values, min_gc)?;

    let style = Into::<ShapeStyle>::into(&RGBColor(0, 200, 0).mix(0.2)).filled();
    plot_collection::mark_subsequences(&mut oe_chart, "O/E", islands, OE_MIN, OE_MAX, style)?;
    plot_collection::mark_subsequences(&mut gc_chart, "GC", islands, GC_MIN, GC_MAX, style)?;

    Ok(())
}

fn process_find_cpgs_cmd(
    fna: &FnaFile,
    args: &FindCpgsArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(charts_dir) = &args.plot {
        let _ = std::fs::create_dir(charts_dir);
    }

    fna.records
        .iter()
        .enumerate()
        .map(|(record_idx, record)| {
            let islands = record
                .content
                .find_cpg_islands(args.min_len, args.step, args.min_gc, args.min_oe);

            let metrics = record.content.get_window_gc_oe_metrics(args.min_len, args.step);

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

                if let Some(charts_dir) = &args.plot {
                    let chart_path = charts_dir.join(format!("{record_idx}.png"));
                    let root =
                        BitMapBackend::new(&chart_path, (1200, 1000)).into_drawing_area();

                    root.fill(&WHITE)?;

                    plot_cpg_islands(
                        &root,
                        &format!("{} #{record_idx}", record.header),
                        args.min_gc,
                        args.min_oe,
                        (&positions, &gc_values, &oe_values),
                        &islands,
                    )?;

                    root.present()?;
                }

                Ok(())
            },
        )
}

fn process_find_orfs_cmd(
    fna: &FnaFile,
    args: &FindOrfsArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    fna.records.iter().enumerate().try_for_each(
        |(record_idx, record)| -> Result<(), Box<dyn std::error::Error>> {
            println!("#{record_idx}: {}.", record.header);

            let orfs = record.content.find_orfs(args.min_len);

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
    )
}

fn process_find_motifs_cmd(
    fna: &FnaFile,
    args: &FindMotifsArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    fna.records.iter().enumerate().try_for_each(
        |(record_id, record)| -> Result<(), Box<dyn std::error::Error>> {
            println!("#{record_id}: {}.", record.header);

            let motif = Sequence::from_str(&args.sub_str)?;

            record
                .content
                .find_motifs(&motif.0)
                .into_iter()
                .for_each(|pos| print!("{pos} "));

            println!();

            Ok(())
        },
    )
}

pub fn process_find_cmd(
    fna: &FnaFile,
    sub_cmd: &FindSubcommand,
) -> Result<(), Box<dyn std::error::Error>> {
    match sub_cmd {
        FindSubcommand::CpgIslands(args) => process_find_cpgs_cmd(fna, args),
        FindSubcommand::Orfs(args) => process_find_orfs_cmd(fna, args),
        FindSubcommand::Motifs(args) => process_find_motifs_cmd(fna, args),
    }
}
