use fna::{FnaFile, FnaRecord};
use gff::{GffFeatureType, GffFile};
use nucleotides::{
    consts,
    sequence::{Sequence, SubSequence},
};

use clap::{Args, Subcommand};
use plot_collection::PlotCpgIslandsOptions;
use plotters::prelude::*;

use std::{path::PathBuf, str::FromStr};

#[derive(Subcommand, Debug)]
pub enum FindSubcommand {
    CpgIslands(FindCpgsArgs),
    Orfs(FindOrfsArgs),
    Motifs(FindMotifsArgs),
}

#[derive(Args, Debug)]
pub struct FindCpgsArgs {
    /// Minimum length of CpG-island
    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_LEN)]
    min_len: usize,

    /// Step of searching CpG-island
    #[arg(long, default_value_t = consts::DEFAULT_CPG_WIN_STEP)]
    step: usize,

    /// Minimum GC% of CpG-island
    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_GC)]
    min_gc: f64,

    /// Minimum O/E of CpG-island
    #[arg(long, default_value_t = consts::DEFAULT_CPG_MIN_OE)]
    min_oe: f64,

    /// Annotate CpG-islands with records of given GFF file (--gff <FILE>)
    #[arg(long)]
    annotate: bool,

    /// Show only main attributes
    #[arg(long)]
    brief: bool,

    /// Draw GC and O/E charts with CpG-islands
    #[arg(long)]
    plot: Option<PathBuf>,
}

#[derive(Args, Debug)]
pub struct FindOrfsArgs {
    /// Minimum ORF length
    #[arg(long, default_value_t = 300)]
    min_len: usize,

    /// Remove ORFs placed in bigger ORFs
    #[arg(long)]
    best_only: bool,

    /// Annotate ORFs with records of given GFF file (--gff <FILE>)
    #[arg(long)]
    annotate: bool,

    /// Show only main attributes
    #[arg(long)]
    brief: bool,
}

#[derive(Args, Debug)]
pub struct FindMotifsArgs {
    #[arg(long)]
    sub_str: String,
}

fn notify_sequences(
    record: &FnaRecord,
    gff: &GffFile,
    subseqs: &[SubSequence],
) -> Vec<(SubSequence, usize, usize)> {
    let features_subseqs = gff
        .features
        .iter()
        .enumerate()
        .filter(|(_, feature)| record.header.starts_with(&feature.seq_id))
        .filter_map(|(feature_id, feature)| {
            let sub_seq =
                SubSequence::new_in(&record.content.0, feature.start - 1, feature.end - 1);

            if sub_seq.is_none() {
                eprintln!("Failed to construct subsequence from feature #{feature_id}!");
            }

            sub_seq
        })
        .collect::<Vec<SubSequence>>();

    record
        .content
        .find_subsequence_intersections(subseqs, &features_subseqs)
}

fn show_features(
    sub_seq_id: usize,
    intersections: &[(SubSequence, usize, usize)],
    gff: &GffFile,
    brief: bool,
) {
    println!("Features:");

    for (_, _, feature_id) in intersections.iter().filter(|(_, id, _)| sub_seq_id == *id) {
        let feature = &gff.features[*feature_id];
        print!(" - {}: ", feature.feature_type);

        if brief {
            match &feature.feature_type {
                GffFeatureType::Gene => {
                    feature
                        .attributes
                        .get("gene")
                        .iter()
                        .for_each(|gene_name| print!("{gene_name} "));
                    feature
                        .attributes
                        .get("locus_tag")
                        .iter()
                        .for_each(|tag| print!("({tag})"));
                }
                GffFeatureType::Cds => {
                    feature
                        .attributes
                        .get("product")
                        .iter()
                        .for_each(|product| print!("{product} "));
                    feature
                        .attributes
                        .get("protein_id")
                        .iter()
                        .for_each(|id| print!("({id})"));
                }
                _ => print!("{}-{}", feature.start - 1, feature.end - 1),
            }
            println!();
        } else {
            println!();
            for (attr_name, attr_value) in feature.attributes.iter() {
                println!("    {attr_name}={attr_value}");
            }
        }
    }
}

fn process_find_cpgs_cmd(
    fna: &FnaFile,
    gff: &Option<GffFile>,
    args: &FindCpgsArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(charts_dir) = &args.plot {
        let _ = std::fs::create_dir(charts_dir);
    }

    if args.annotate && gff.is_none() {
        return Err("--annotate requires specified GFF (--gff <FILE>)".into());
    }

    fna.records
        .iter()
        .enumerate()
        .map(|(record_idx, record)| {
            let islands =
                record
                    .content
                    .find_cpg_islands(args.min_len, args.step, args.min_gc, args.min_oe);

            let metrics = record
                .content
                .get_window_gc_oe_metrics(args.min_len, args.step);

            (record_idx, record, islands, metrics)
        })
        .try_for_each(
            |(record_idx, record, islands, metrics)| -> Result<(), Box<dyn std::error::Error>> {
                println!("#{record_idx}: {}.", record.header);

                let intersections = if args.annotate {
                    Some(notify_sequences(record, gff.as_ref().unwrap(), &islands))
                } else {
                    None
                };

                for (island_id, island) in islands.iter().enumerate() {
                    println!(
                        "\nIsland #{island_id}:\nposition: {}-{} (total len: {})\nGC: {}%",
                        island.start,
                        island.end,
                        island.seq.length(),
                        island.seq.get_gc_content()
                    );

                    if let Some(intersections) = &intersections {
                        show_features(island_id, intersections, gff.as_ref().unwrap(), args.brief);
                    }
                }

                if let Some(charts_dir) = &args.plot {
                    let chart_path = charts_dir.join(format!("{record_idx}.png"));
                    let root = BitMapBackend::new(&chart_path, (1200, 1000)).into_drawing_area();

                    root.fill(&WHITE)?;

                    let opts = PlotCpgIslandsOptions {
                        show_gc: true,
                        show_oe: true,
                        show_min_gc: true,
                        show_min_oe: true,
                        show_cpg_islands: true,
                        show_orfs: false,
                        min_gc: args.min_gc,
                        min_oe: args.min_oe,
                    };

                    let title = format!("{} #{record_idx}", record.header);

                    plot_collection::plot_cpg_islands(&root, &title, &metrics, &islands, &opts)?;

                    root.present()?;
                }

                Ok(())
            },
        )
}

fn process_find_orfs_cmd(
    fna: &FnaFile,
    gff: &Option<GffFile>,
    args: &FindOrfsArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    if args.annotate && gff.is_none() {
        return Err("--annotate requires specified GFF (--gff <FILE>)".into());
    }

    fna.records.iter().enumerate().try_for_each(
        |(record_idx, record)| -> Result<(), Box<dyn std::error::Error>> {
            println!("#{record_idx}: {}.", record.header);

            let orfs = record.content.find_orfs(args.min_len);
            let orfs = Sequence::best_orfs(orfs);

            let intersections = if args.annotate {
                Some(notify_sequences(record, gff.as_ref().unwrap(), &orfs))
            } else {
                None
            };

            for (orf_id, orf) in orfs.iter().enumerate() {
                println!(
                    "\nORF #{orf_id}:\nposition: {}-{} (total len: {})\nGC: {}%",
                    orf.start,
                    orf.end,
                    orf.seq.length(),
                    orf.seq.get_gc_content()
                );

                if let Some(intersections) = &intersections {
                    show_features(orf_id, intersections, gff.as_ref().unwrap(), args.brief);
                }
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
    gff: &Option<GffFile>,
    sub_cmd: &FindSubcommand,
) -> Result<(), Box<dyn std::error::Error>> {
    match sub_cmd {
        FindSubcommand::CpgIslands(args) => process_find_cpgs_cmd(fna, gff, args),
        FindSubcommand::Orfs(args) => process_find_orfs_cmd(fna, gff, args),
        FindSubcommand::Motifs(args) => process_find_motifs_cmd(fna, args),
    }
}
