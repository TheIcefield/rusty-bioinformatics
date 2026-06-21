use nucleotides::sequence::SubSequence;
use plotters::{
    backend::BitMapBackend,
    chart::ChartContext,
    coord::{cartesian::Cartesian2d, types::RangedCoordf64},
    element::{PathElement, Rectangle},
    series::LineSeries,
    style::{
        RGBColor, ShapeStyle,
        full_palette::{BLUE, RED},
    },
};

pub fn plot_gc_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    gc_values: &[f64],
    min_gc: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    const GC_VALUES_COLOR: RGBColor = BLUE;
    const MIN_GC_COLOR: RGBColor = RED;

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

    chart
        .draw_series(LineSeries::new(
            (x_min as i32..=x_max as i32).map(|x| (x as f64, min_gc)),
            &MIN_GC_COLOR,
        ))?
        .label("Min GC %")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], MIN_GC_COLOR));

    Ok(())
}

pub fn plot_oe_values(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    positions: &[usize],
    oe_values: &[f64],
    min_oe: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    const OE_VALUES_COLOR: RGBColor = BLUE;
    const MIN_OE_COLOR: RGBColor = RED;

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

    chart
        .draw_series(LineSeries::new(
            (x_min as i32..=x_max as i32).map(|x| (x as f64, min_oe)),
            &MIN_OE_COLOR,
        ))?
        .label("Min O/E %")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], MIN_OE_COLOR));

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
