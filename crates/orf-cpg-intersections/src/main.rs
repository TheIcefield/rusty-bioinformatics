use clap::Parser;

use fna::FnaFile;
use nucleotides::sequence::{CpgIsland, Orf, Sequence};
use plotters::{coord::Shift, prelude::*};

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

fn plot_cpg_orf_intersections(
    root: &DrawingArea<BitMapBackend, Shift>,
    title: &str,
    min_gc: f64,
    min_oe: f64,
    metrics: (&[usize], &[f64], &[f64]),
    islands: &[CpgIsland],
    orfs: &[Orf],
) -> Result<(), Box<dyn std::error::Error>> {
    const FONT_NAME: &str = "sans-serif";
    const TITLE_FONT_SIZE: f64 = 30f64;
    const CPG_ISLAND_LEGEND: &str = "CpG-islands";
    const ORF_LEGEND: &str = "ORF";
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

    let orf_style = Into::<ShapeStyle>::into(&RGBColor(150, 100, 0).mix(0.2)).filled();
    let island_style = Into::<ShapeStyle>::into(&RGBColor(0, 200, 0).mix(0.3)).filled();

    // Draw ORF
    plot_collection::mark_subsequences(&mut gc_chart, ORF_LEGEND, orfs, GC_MIN, GC_MAX, orf_style)?;
    plot_collection::mark_subsequences(&mut oe_chart, ORF_LEGEND, orfs, OE_MIN, OE_MAX, orf_style)?;

    // Draw CpG-islands
    plot_collection::mark_subsequences(
        &mut gc_chart,
        CPG_ISLAND_LEGEND,
        islands,
        GC_MIN,
        GC_MAX,
        island_style,
    )?;

    plot_collection::mark_subsequences(
        &mut oe_chart,
        CPG_ISLAND_LEGEND,
        islands,
        OE_MIN,
        OE_MAX,
        island_style,
    )?;

    // Draw GC and O/E values
    plot_collection::plot_gc_values(&mut gc_chart, positions, gc_values, min_gc)?;
    plot_collection::plot_oe_values(&mut oe_chart, positions, oe_values, min_oe)?;

    oe_chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    gc_chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
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
        .filter_map(|(record_id, record)| {
            let islands = record.content.find_cpg_islands(
                cli.cgp_min_len,
                cli.cpg_min_gc,
                cli.cpg_min_oe,
                cli.cgp_step,
            );

            let Ok(orfs) = record.content.find_orfs(cli.orf_min_len) else {
                return None;
            };

            let orfs = if cli.best_orf_only { Sequence::best_orfs(orfs) } else { orfs };

            let intersections = record.content.find_subsequence_intersections(&islands, &orfs);

            let metrics = record.content.get_window_gc_oe_metrics(std::cmp::min(cli.cgp_min_len, cli.orf_min_len), cli.cgp_step);

            Some((record_id, record, islands, orfs, intersections, metrics))
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

                    plot_cpg_orf_intersections(
                        &root,
                        &format!("{} #{record_id}", record.header),
                        cli.cpg_min_gc,
                        cli.cpg_min_oe,
                        (&metrics.0, &metrics.1, &metrics.2),
                        &islands,
                        &orfs
                    )?;

                    root.present()?;
                }

                Ok(())
            },
        )?;

    Ok(())
}
