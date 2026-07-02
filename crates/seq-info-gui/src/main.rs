use std::path::Path;

use fna::FnaFile;
use iced::{
    Alignment, Color, Element, Length, Pixels,
    widget::{button, column, container, pick_list, row, text, toggler},
};
use iced_plot::{LineStyle, PlotUiMessage, PlotWidget, PlotWidgetBuilder, Series};
use nucleotides::{
    consts,
    sequence::{CpgIsland, GcMetrics, Sequence},
};
use plot_collection::PlotCpgIslandsOptions;
use plotters::{backend::BitMapBackend, drawing::IntoDrawingArea, style::full_palette::WHITE};

const GC_PLOT_ID: usize = 1;
const OE_PLOT_ID: usize = 2;

const TITLE_FONT_SIZE: f32 = 25.;
const GC_PLOT_TITLE: &str = "GC percentage";
const OE_PLOT_TITLE: &str = "Observed / Expected CG motifs";

fn main() -> iced::Result {
    iced::application(AppState::new, AppState::update, AppState::view).run()
}

#[derive(Debug, Clone)]
enum Message {
    LoadFile,
    PlotUpdate { msg: PlotUiMessage, plot_id: usize },
    SavePlots,
    StrandSelected(usize),
    ViewOptionsChanged(PlotCpgIslandsOptions),
}

#[derive(Default)]
struct AppState {
    current_strand: Option<usize>,
    strands: Vec<SequenceInfo>,
    view: PlotCpgIslandsOptions,
}

#[derive(Default)]
struct SequenceInfo {
    header: String,
    gc_content_metrics: Option<GcMetrics>,
    cpg_islands: Option<Vec<CpgIsland>>,
    orfs: Option<Vec<CpgIsland>>,
    sequence: Sequence,
    gc_plot: PlotWidget,
    oe_plot: PlotWidget,
}

impl SequenceInfo {
    fn update(&mut self, options: &PlotCpgIslandsOptions) {
        if (options.show_gc || options.show_oe) && self.gc_content_metrics.is_none() {
            self.collect_metrics();
        }
        if options.show_cpg_islands && self.cpg_islands.is_none() {
            self.collect_cpg_islands();
        }
        if options.show_orfs && self.orfs.is_none() {
            self.collect_orfs();
        }

        self.build_gc_plot(options);
        self.build_oe_plot(options);
    }

    fn build_gc_plot(&mut self, options: &PlotCpgIslandsOptions) {
        let mut temp_plot = PlotWidgetBuilder::new()
            .with_x_label("Nucleotide")
            .with_y_label("GC%")
            .with_x_tick_labels(true)
            .with_y_tick_labels(true)
            .with_cursor_overlay(true)
            .with_crosshairs(true);

        let metrics = self.gc_content_metrics.as_ref().unwrap();

        if options.show_gc {
            let gc_value = metrics
                .positions
                .iter()
                .zip(metrics.gc_values.iter())
                .map(|(pos, gc)| [*pos as f64, *gc])
                .collect();

            temp_plot = temp_plot.add_series(
                Series::circles(gc_value, 1.0)
                    .line_solid()
                    .with_label("GC %"),
            );

            if options.show_min_gc {
                temp_plot = temp_plot.add_series(
                    Series::line_only(
                        metrics
                            .positions
                            .iter()
                            .map(|x| [*x as f64, consts::DEFAULT_CPG_MIN_GC])
                            .collect(),
                        LineStyle::solid(),
                    )
                    .line_solid()
                    .with_label("Min GC%")
                    .with_color(Color::from_rgb(1.0, 0.0, 0.0)),
                )
            }
        }

        self.gc_plot = temp_plot.build().unwrap();
    }

    fn build_oe_plot(&mut self, options: &PlotCpgIslandsOptions) {
        let mut temp_plot = PlotWidgetBuilder::new()
            .with_x_label("Nucleotide")
            .with_y_label("O/E")
            .with_x_tick_labels(true)
            .with_y_tick_labels(true)
            .with_cursor_overlay(true)
            .with_crosshairs(true);

        let metrics = self.gc_content_metrics.as_ref().unwrap();

        if options.show_oe {
            let oe_value = metrics
                .positions
                .iter()
                .zip(metrics.oe_values.iter())
                .map(|(pos, oe)| [*pos as f64, *oe])
                .collect();

            temp_plot = temp_plot.add_series(
                Series::circles(oe_value, 1.0)
                    .line_solid()
                    .with_label("O/E"),
            );

            if options.show_min_oe {
                temp_plot = temp_plot.add_series(
                    Series::line_only(
                        metrics
                            .positions
                            .iter()
                            .map(|x| [*x as f64, consts::DEFAULT_CPG_MIN_OE])
                            .collect(),
                        LineStyle::solid(),
                    )
                    .line_solid()
                    .with_label("Min O/E")
                    .with_color(Color::from_rgb(1.0, 0.0, 0.0)),
                )
            }
        }

        self.oe_plot = temp_plot.build().unwrap();
    }

    fn collect_metrics(&mut self) {
        let metrics = self
            .sequence
            .get_window_gc_oe_metrics(consts::DEFAULT_CPG_MIN_LEN, consts::DEFAULT_CPG_WIN_STEP);
        self.gc_content_metrics = Some(metrics);
    }

    fn collect_cpg_islands(&mut self) {
        let islands = self.sequence.find_cpg_islands(
            consts::DEFAULT_CPG_MIN_LEN,
            consts::DEFAULT_CPG_WIN_STEP,
            consts::DEFAULT_CPG_MIN_GC,
            consts::DEFAULT_CPG_MIN_OE,
        );
        self.cpg_islands = Some(islands);
    }

    fn collect_orfs(&mut self) {
        let orfs = self.sequence.find_orfs(consts::DEFAULT_ORF_MIN_LEN);
        self.orfs = Some(orfs);
    }
}

impl AppState {
    fn new() -> Self {
        Self {
            current_strand: None,
            strands: Vec::new(),
            view: PlotCpgIslandsOptions {
                show_gc: true,
                show_oe: true,
                min_gc: consts::DEFAULT_CPG_MIN_GC,
                min_oe: consts::DEFAULT_CPG_MIN_OE,
                ..Default::default()
            },
        }
    }

    fn update(&mut self, message: Message) {
        match message {
            Message::LoadFile => {
                if let Some(path) = rfd::FileDialog::new().pick_file() {
                    self.load_fna(&path);
                }
            }
            Message::PlotUpdate { msg, plot_id } => {
                if let Some(strand) = self.get_choosen_strand_mut() {
                    match plot_id {
                        GC_PLOT_ID => strand.gc_plot.update(msg),
                        OE_PLOT_ID => strand.oe_plot.update(msg),
                        _ => {}
                    }
                }
            }
            Message::SavePlots => {
                if let Some(path) = rfd::FileDialog::new().save_file()
                    && let Err(err) = self.save_plots(&path)
                {
                    eprintln!("{err}");
                }
            }
            Message::StrandSelected(index) => {
                self.current_strand = Some(index);
                let options = self.view;

                if let Some(strand) = self.get_choosen_strand_mut() {
                    strand.update(&options);
                }
            }
            Message::ViewOptionsChanged(options) => {
                self.view = options;

                if let Some(strand) = self.get_choosen_strand_mut() {
                    strand.update(&options);
                }
            }
        }
    }

    fn load_fna(&mut self, path: &Path) {
        let fna = FnaFile::open(path).unwrap();

        for record in fna.records {
            let mut seq = SequenceInfo {
                header: record.header,
                sequence: record.content,
                gc_content_metrics: None,
                cpg_islands: None,
                orfs: None,
                gc_plot: PlotWidgetBuilder::new().build().unwrap(),
                oe_plot: PlotWidgetBuilder::new().build().unwrap(),
            };
            seq.update(&self.view);
            self.strands.push(seq);
        }

        // Choose last loaded
        self.current_strand = Some(self.strands.len() - 1);
    }

    fn save_plots(&mut self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let opts = self.view;

        let Some(strand) = self.get_choosen_strand_mut() else {
            return Err("No choosen strand".into());
        };

        strand.update(&opts);

        let metrics = strand.gc_content_metrics.as_ref().unwrap();

        let islands = if opts.show_cpg_islands {
            Some(strand.cpg_islands.as_ref().unwrap().as_ref())
        } else {
            None
        };

        let orfs = if opts.show_orfs {
            Some(strand.orfs.as_ref().unwrap().as_ref())
        } else {
            None
        };

        let root = BitMapBackend::new(path, (1200, 1000)).into_drawing_area();
        root.fill(&WHITE)?;

        plot_collection::plot_cpg_and_orf(&root, &strand.header, metrics, &islands, &orfs, &opts)?;

        root.present()?;

        Ok(())
    }

    fn view_plot_panel<'a>(title: &str, plot: Element<'a, Message>) -> Element<'a, Message> {
        column![
            text!("{title}"),
            container(plot).width(Length::Fill).height(Length::Fill)
        ]
        .spacing(8)
        .align_x(Alignment::Center)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
    }

    fn view_checkbox_option<'a, F>(title: &str, value: bool, on_toggle: F) -> Element<'a, Message>
    where
        F: 'a + Fn(bool) -> Message,
    {
        column![text!("{title}"), toggler(value).on_toggle(on_toggle)]
            .align_x(Alignment::Center)
            .into()
    }

    fn view_checkbox_options<'a>(current_view: &'a PlotCpgIslandsOptions) -> Element<'a, Message> {
        row![
            Self::view_checkbox_option("Show GC", current_view.show_gc, |options| {
                let mut modified = *current_view;
                modified.show_gc = options;
                Message::ViewOptionsChanged(modified)
            }),
            Self::view_checkbox_option("Show min GC", current_view.show_min_gc, |options| {
                let mut modified = *current_view;
                modified.show_min_gc = options;
                Message::ViewOptionsChanged(modified)
            }),
            Self::view_checkbox_option("Show O/E", current_view.show_oe, |options| {
                let mut modified = *current_view;
                modified.show_oe = options;
                Message::ViewOptionsChanged(modified)
            }),
            Self::view_checkbox_option("Show min O/E", current_view.show_min_oe, |options| {
                let mut modified = *current_view;
                modified.show_min_oe = options;
                Message::ViewOptionsChanged(modified)
            }),
            Self::view_checkbox_option(
                "Show CpG-islands",
                current_view.show_cpg_islands,
                |options| {
                    let mut modified = *current_view;
                    modified.show_cpg_islands = options;
                    Message::ViewOptionsChanged(modified)
                }
            ),
            Self::view_checkbox_option("Show ORFs", current_view.show_orfs, |options| {
                let mut modified = *current_view;
                modified.show_orfs = options;
                Message::ViewOptionsChanged(modified)
            })
        ]
        .spacing(10)
        .into()
    }

    fn view(&self) -> Element<'_, Message> {
        let mut plots = column![];

        if let Some(strand) = self.get_choosen_strand() {
            plots = plots.push(text(strand.header.as_str()).size(Pixels::from(TITLE_FONT_SIZE)));

            if self.view.show_gc {
                plots = plots.push(Self::view_plot_panel(
                    GC_PLOT_TITLE,
                    strand.gc_plot.view().map(|msg| Message::PlotUpdate {
                        msg,
                        plot_id: GC_PLOT_ID,
                    }),
                ));
            }

            if self.view.show_oe {
                plots = plots.push(Self::view_plot_panel(
                    OE_PLOT_TITLE,
                    strand.oe_plot.view().map(|msg| Message::PlotUpdate {
                        msg,
                        plot_id: OE_PLOT_ID,
                    }),
                ));
            }
        }

        column![
            row![
                button("Load FNA").on_press(Message::LoadFile),
                button("Save Plots").on_press(Message::SavePlots),
                pick_list(
                    (0..self.strands.len()).collect::<Vec<usize>>(),
                    self.current_strand,
                    Message::StrandSelected
                )
                .placeholder("Choose strand"),
                Self::view_checkbox_options(&self.view),
            ]
            .spacing(10)
            .align_y(Alignment::Center),
            plots.align_x(Alignment::Center).spacing(5)
        ]
        .spacing(10)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
    }

    fn get_choosen_strand(&self) -> Option<&SequenceInfo> {
        let strand_id = self.current_strand?;
        let strand = self.strands.get(strand_id)?;

        Some(strand)
    }

    fn get_choosen_strand_mut(&mut self) -> Option<&mut SequenceInfo> {
        let strand_id = self.current_strand?;
        let strand = self.strands.get_mut(strand_id)?;

        Some(strand)
    }
}
