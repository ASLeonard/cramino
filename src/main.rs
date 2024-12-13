use clap::Parser;
use std::path::PathBuf;

pub mod calculations;
pub mod extract_from_bam;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Tool to extract QC metrics from cram or bam", long_about = None)]
pub struct Cli {
    /// cram or bam file to check
    #[clap(value_parser, default_value = "-")]
    input: String,

    /// Number of parallel decompression threads to use
    #[clap(short, long, value_parser, default_value_t = 4)]
    threads: usize,

    /// reference for decompressing cram
    #[clap(long, value_parser)]
    reference: Option<String>,

}

pub fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() -> Result<(), rust_htslib::errors::Error> {
   let args = Cli::parse(); 
   let (metrics, header) = extract_from_bam::extract(&args);

   let lengths = metrics.lengths.as_ref().unwrap();
   let qualities = metrics.qualities.as_ref().unwrap();

   let rg = &header.to_hashmap()["RG"];
   if rg.len() > 1 {
     panic!("Multiple read groups present in the file!");
   }
   
   let data_yield = lengths.iter().sum::<u64>();
   println!("{},{},{},{},{:.3}",rg[0]["PU"],lengths.len(), data_yield, calculations::get_n(lengths, data_yield, 0.50), calculations::mean_accuracy(qualities)) ;


   Ok(())
}
