use nucleotides::sequence::CpgIsland;

use clap::Subcommand;
use plotters::{coord::Shift, prelude::*};

use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub enum FindSubcommand {
    CpgIslands {
        #[arg(long, default_value_t = 200)]
        min_len: usize,

        #[arg(long, default_value_t = 50.0)]
        min_gc: f64,

        #[arg(long, default_value_t = 0.6)]
        min_oe: f64,

        #[arg(long, default_value_t = 10)]
        step: usize,

        #[arg(long)]
        plot: Option<PathBuf>,
    },

    Orf {
        #[arg(long, default_value_t = 300)]
        min_len: usize,
    },

    Motif {
        #[arg(long)]
        sub_str: String,
    },
}

pub fn plot_cpg_islands(
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
