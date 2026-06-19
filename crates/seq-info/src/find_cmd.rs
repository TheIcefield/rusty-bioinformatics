use nucleotides::sequence::CpgIsland;

use clap::Subcommand;
use plotters::{
    coord::{Shift, types::RangedCoordf64},
    prelude::*,
};

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
}

fn plot_gc_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    gc_values: &[f64],
    min_gc: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let x_min = *positions.first().unwrap_or(&0) as f64;
    let x_max = *positions.last().unwrap_or(&0) as f64;

    chart.configure_mesh().x_labels(20).y_labels(10).draw()?;

    chart
        .draw_series(LineSeries::new(
            positions
                .iter()
                .zip(gc_values.iter())
                .map(|(&x, &y)| (x as f64, y)),
            &BLUE,
        ))?
        .label("GC percentage");

    chart
        .draw_series(LineSeries::new(
            (x_min as i32..=x_max as i32).map(|x| (x as f64, min_gc)),
            &RED,
        ))?
        .label("Min GC%");

    Ok(())
}

fn plot_oe_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    oe_values: &[f64],
    min_oe: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let x_min = *positions.first().unwrap_or(&0) as f64;
    let x_max = *positions.last().unwrap_or(&0) as f64;

    chart.configure_mesh().x_labels(20).y_labels(10).draw()?;

    chart.draw_series(LineSeries::new(
        positions
            .iter()
            .zip(oe_values.iter())
            .map(|(&x, &y)| (x as f64, y)),
        &BLUE,
    ))?;

    chart.draw_series(LineSeries::new(
        (x_min as i32..=x_max as i32).map(|x| (x as f64, min_oe)),
        &RED,
    ))?;

    Ok(())
}

fn mark_cpg_islands(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    islands: &[CpgIsland],
    min: f64,
    max: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    for island in islands {
        let style = Into::<ShapeStyle>::into(&RGBColor(0, 200, 0).mix(0.2)).filled();
        let rect = Rectangle::new(
            [(island.start as f64, min), (island.end as f64, max)],
            style,
        );
        let rect_series = vec![rect];

        chart.draw_series(rect_series)?;
    }

    Ok(())
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

    let mut upper_chart = ChartBuilder::on(&upper)
        .caption("Observed/Expected CpG", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .margin_top(30)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, OE_MIN..OE_MAX)?;

    let mut lower_chart = ChartBuilder::on(&lower)
        .caption("GC content (%)", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, GC_MIN..GC_MAX)?;

    plot_gc_values(&mut upper_chart, positions, oe_values, min_oe)?;
    plot_oe_values(&mut lower_chart, positions, gc_values, min_gc)?;

    mark_cpg_islands(&mut upper_chart, islands, OE_MIN, OE_MAX)?;
    mark_cpg_islands(&mut lower_chart, islands, GC_MIN, GC_MAX)?;

    root.present()?;

    Ok(())
}
