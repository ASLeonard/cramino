use clap::Parser;
use extract_from_bam::Data;
use std::path::PathBuf;

pub mod calculations;
pub mod extract_from_bam;
//pub mod feather;
//pub mod file_info;
//pub mod histograms;
//pub mod karyotype;
//pub mod phased;
//pub mod splicing;
pub mod utils;

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
  
   for (key, value) in header.to_hashmap().into_iter() {
      println!("{} / {:?}", key, value);
   }

   metrics_from_bam(metrics)?;
   Ok(())
}

fn metrics_from_bam(
    metrics: Data,
) -> Result<(), rust_htslib::errors::Error> {

    generate_main_output(
        metrics.lengths.as_ref().unwrap(),
        metrics.qualities.as_ref().unwrap(),
        metrics.all_counts,
    );

    Ok(())
}

fn generate_main_output(
    lengths: &Vec<u128>,
    qualities: &Vec<f64>,
    all_reads: usize,
) {
    let num_reads = lengths.len();
    let data_yield: u128 = lengths.iter().sum::<u128>();
    println!("Number of alignments\t{num_reads}");
    println!(
        "% from total reads\t{:.2}",
        (num_reads as f64) / (all_reads as f64) * 100.0
    );
    println!("Yield [Gb]\t{:.2}", data_yield as f64 / 1e9);
    let data_yield_long = lengths.iter().filter(|l| l > &&25000).sum::<u128>();
    println!("Yield [Gb] (>25kb)\t{:.2}", data_yield_long as f64 / 1e9);


    println!("N50\t{}", calculations::get_n(lengths, data_yield, 0.50));
    println!("QV \t{:.2}", qualities.iter().sum::<f64>() / qualities.len() as f64);   
 
}
