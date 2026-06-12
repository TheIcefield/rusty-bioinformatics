use std::{io::Read, path::Path};

use fna::FnaFile;

fn help() -> Result<(), String> {
    println!(
        "Help:\
       \nnucl-gc-content <path/to/file.fna>"
    );

    Ok(())
}

fn get_gc_content(path: &Path) -> Result<(), String> {
    let mut file = std::fs::File::open(path).map_err(|err| err.to_string())?;

    let mut raw_data = String::new();
    file.read_to_string(&mut raw_data)
        .map_err(|err| err.to_string())?;

    let fna = FnaFile::read(&raw_data)?;

    for (idx, record) in fna.records.iter().enumerate() {
        println!(
            "#{idx}: {}.\n    GC: {}%",
            record.header,
            record.content.get_gc_content()
        )
    }

    Ok(())
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();

    if args.len() != 1 {
        return Err("Expected one argument".to_string());
    }

    match &args[0][..] {
        "--help" | "-h" => help(),
        path => get_gc_content(Path::new(path)),
    }
}
