use clap::Subcommand;

use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub enum CompareSubcommand {
    HammingDistance {
        #[arg(long)]
        second: PathBuf,
    },
}
