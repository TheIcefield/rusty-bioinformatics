use nucleotides::sequence::{CpgIsland, GcMetrics, Orf, SubSequence};
use plotters::{
    backend::BitMapBackend,
    chart::{ChartBuilder, ChartContext},
    coord::{Shift, cartesian::Cartesian2d, types::RangedCoordf64},
    drawing::DrawingArea,
    element::{PathElement, Rectangle},
    series::LineSeries,
    style::{
        Color, RGBColor, ShapeStyle,
        full_palette::{BLACK, BLUE, RED, WHITE},
    },
};

const OE_VALUES_COLOR: RGBColor = BLUE;
const GC_VALUES_COLOR: RGBColor = BLUE;

const MIN_GC_COLOR: RGBColor = RED;
const MIN_OE_COLOR: RGBColor = RED;

#[derive(Default, Debug, Clone, Copy)]
pub struct PlotCpgIslandsOptions {
    pub show_gc: bool,
    pub show_oe: bool,
    pub show_min_gc: bool,
    pub show_min_oe: bool,
    pub show_cpg_islands: bool,
    pub show_orfs: bool,
    pub min_gc: f64,
    pub min_oe: f64,
}

pub fn plot_gc_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    gc_values: &[f64],
    options: &PlotCpgIslandsOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    if options.show_gc {
        let x_min = *positions.first().unwrap_or(&0) as f64;
        let x_max = *positions.last().unwrap_or(&0) as f64;

        chart.configure_mesh().x_labels(20).y_labels(10).draw()?;

        chart
            .draw_series(LineSeries::new(
                positions
                    .iter()
                    .zip(gc_values.iter())
                    .map(|(&x, &y)| (x as f64, y)),
                &GC_VALUES_COLOR,
            ))?
            .label("GC %")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GC_VALUES_COLOR));

        if options.show_min_gc {
            chart
                .draw_series(LineSeries::new(
                    (x_min as i32..=x_max as i32).map(|x| (x as f64, options.min_gc)),
                    &MIN_GC_COLOR,
                ))?
                .label("Min GC %")
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], MIN_GC_COLOR));
        }
    }

    Ok(())
}

pub fn plot_oe_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    oe_values: &[f64],
    options: &PlotCpgIslandsOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    if options.show_oe {
        let x_min = *positions.first().unwrap_or(&0) as f64;
        let x_max = *positions.last().unwrap_or(&0) as f64;

        chart.configure_mesh().x_labels(20).y_labels(10).draw()?;

        chart
            .draw_series(LineSeries::new(
                positions
                    .iter()
                    .zip(oe_values.iter())
                    .map(|(&x, &y)| (x as f64, y)),
                &OE_VALUES_COLOR,
            ))?
            .label("O/E %")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], OE_VALUES_COLOR));

        if options.show_min_oe {
            chart
                .draw_series(LineSeries::new(
                    (x_min as i32..=x_max as i32).map(|x| (x as f64, options.min_oe)),
                    &MIN_OE_COLOR,
                ))?
                .label("Min O/E %")
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], MIN_OE_COLOR));
        }
    }
    Ok(())
}

pub fn mark_subsequences(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    label: &str,
    subsequences: &[SubSequence],
    min: f64,
    max: f64,
    style: ShapeStyle,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut rect_series = Vec::with_capacity(subsequences.len());

    for subseq in subsequences {
        let rect = Rectangle::new(
            [(subseq.start as f64, min), (subseq.end as f64, max)],
            style,
        );
        rect_series.push(rect);
    }

    chart
        .draw_series(rect_series)?
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style));

    Ok(())
}

pub fn plot_cpg_islands(
    root: &DrawingArea<BitMapBackend, Shift>,
    title: &str,
    metrics: &GcMetrics,
    islands: &[CpgIsland],
    options: &PlotCpgIslandsOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    const FONT_NAME: &str = "sans-serif";
    const TITLE_FONT_SIZE: f64 = 30f64;
    const CAPTION_FONT_SIZE: f64 = 20f64;
    const OE_MIN: f64 = 0.0;
    const OE_MAX: f64 = 2.0;
    const GC_MIN: f64 = 0.0;
    const GC_MAX: f64 = 100.0;

    let x_min = *metrics.positions.first().unwrap_or(&0) as f64;
    let x_max = *metrics.positions.last().unwrap_or(&0) as f64;

    root.titled(title, (FONT_NAME, TITLE_FONT_SIZE))?;

    let (upper, lower) = root.split_vertically(500);

    let mut gc_chart = ChartBuilder::on(&upper)
        .caption("GC content (%)", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, GC_MIN..GC_MAX)?;

    let mut oe_chart = ChartBuilder::on(&lower)
        .caption("Observed/Expected CpG", (FONT_NAME, CAPTION_FONT_SIZE))
        .margin(10)
        .margin_top(30)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, OE_MIN..OE_MAX)?;

    plot_gc_values(
        &mut gc_chart,
        &metrics.positions,
        &metrics.gc_values,
        options,
    )?;

    plot_oe_values(
        &mut oe_chart,
        &metrics.positions,
        &metrics.oe_values,
        options,
    )?;

    let style = Into::<ShapeStyle>::into(&RGBColor(0, 200, 0).mix(0.2)).filled();
    mark_subsequences(&mut oe_chart, "O/E", islands, OE_MIN, OE_MAX, style)?;
    mark_subsequences(&mut gc_chart, "GC", islands, GC_MIN, GC_MAX, style)?;

    Ok(())
}

pub fn plot_cpg_and_orf(
    root: &DrawingArea<BitMapBackend, Shift>,
    title: &str,
    metrics: &GcMetrics,
    islands: &Option<&[CpgIsland]>,
    orfs: &Option<&[Orf]>,
    options: &PlotCpgIslandsOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    const FONT_NAME: &str = "sans-serif";
    const CPG_ISLAND_LEGEND: &str = "CpG-islands";
    const ORF_LEGEND: &str = "ORF";
    const TITLE_FONT_SIZE: f64 = 30f64;
    const CAPTION_FONT_SIZE: f64 = 20f64;
    const OE_MIN: f64 = 0.0;
    const OE_MAX: f64 = 2.0;
    const GC_MIN: f64 = 0.0;
    const GC_MAX: f64 = 100.0;

    let x_min = *metrics.positions.first().unwrap_or(&0) as f64;
    let x_max = *metrics.positions.last().unwrap_or(&0) as f64;

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

    // Draw ORF
    if options.show_orfs
        && let Some(orfs) = orfs
    {
        let orf_style = Into::<ShapeStyle>::into(&RGBColor(150, 100, 0).mix(0.2)).filled();

        mark_subsequences(&mut gc_chart, ORF_LEGEND, orfs, GC_MIN, GC_MAX, orf_style)?;
        mark_subsequences(&mut oe_chart, ORF_LEGEND, orfs, OE_MIN, OE_MAX, orf_style)?;
    }

    // Draw CpG-islands
    if options.show_cpg_islands
        && let Some(islands) = islands
    {
        let island_style = Into::<ShapeStyle>::into(&RGBColor(0, 200, 0).mix(0.3)).filled();

        mark_subsequences(
            &mut gc_chart,
            CPG_ISLAND_LEGEND,
            islands,
            GC_MIN,
            GC_MAX,
            island_style,
        )?;

        mark_subsequences(
            &mut oe_chart,
            CPG_ISLAND_LEGEND,
            islands,
            OE_MIN,
            OE_MAX,
            island_style,
        )?;
    }

    // Draw GC and O/E values
    plot_gc_values(
        &mut gc_chart,
        &metrics.positions,
        &metrics.gc_values,
        options,
    )?;

    plot_oe_values(
        &mut oe_chart,
        &metrics.positions,
        &metrics.oe_values,
        options,
    )?;

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
