use clap::Args;
use fna::FnaFile;
use nucleotides::nucleotide::Nucleotide;

#[derive(Args, Debug)]
pub struct InfoSubcommandArgs {
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
}

pub fn process_info_cmd(
    fna: &FnaFile,
    args: &InfoSubcommandArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    for (idx, record) in fna.records.iter().enumerate() {
        println!("#{idx}: {}.", record.header);

        if args.kind {
            println!("    Sequence kind: {}", record.content.get_kind());
        }

        if args.len {
            println!("    Length: {} nucleotides", record.content.length());
        }

        if args.count {
            println!("    A: {}", record.content.count(Nucleotide::A));

            if record.content.is_dna() {
                println!("    T: {}", record.content.count(Nucleotide::T));
            } else {
                println!("    T: {}", record.content.count(Nucleotide::U));
            }

            println!("    G: {}", record.content.count(Nucleotide::G));
            println!("    C: {}", record.content.count(Nucleotide::C));
        }

        if args.gc {
            println!("    GC content: {}%", record.content.get_gc_content());
        }
    }

    Ok(())
}
