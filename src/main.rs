use clap::Parser;
use std::path::{Path, PathBuf};

pub mod calculations;
pub mod extract_from_bam;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Tool to extract QC metrics from cram or bam", long_about = None)]
pub struct Cli {
    /// cram or bam file to check
    #[clap(value_parser, default_value = "-")]
    pub input: String,

    /// Number of parallel decompression threads to use
    #[clap(short, long, value_parser, default_value_t = 4)]
    pub threads: usize,

    /// reference for decompressing cram
    #[clap(long, value_parser)]
    pub reference: Option<String>,
}

pub fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn extract_pu(header: &rust_htslib::bam::Header, input_path: &str) -> String {
    let header_map = header.to_hashmap();
    if let Some(rg_records) = header_map.get("RG") {
        let mut pus = Vec::new();
        for record in rg_records {
            if let Some(pu) = record.get("PU").map(String::as_str).filter(|s| !s.is_empty()) {
                if !pus.contains(&pu) {
                    pus.push(pu);
                }
            }
        }
        if !pus.is_empty() {
            return pus.join(";");
        }

        let mut ids = Vec::new();
        for record in rg_records {
            if let Some(id) = record.get("ID").map(String::as_str).filter(|s| !s.is_empty()) {
                if !ids.contains(&id) {
                    ids.push(id);
                }
            }
        }
        if !ids.is_empty() {
            return ids.join(";");
        }
    }

    if input_path == "-" {
        "stdin".to_string()
    } else {
        Path::new(input_path)
            .file_name()
            .and_then(|name| name.to_str())
            .filter(|name| !name.is_empty())
            .unwrap_or("NA")
            .to_string()
    }
}

fn main() -> Result<(), rust_htslib::errors::Error> {
    let args = Cli::parse();
    let (metrics, header) = extract_from_bam::extract(&args);

    let pu = extract_pu(&header, &args.input);

    println!(
        "{},{},{},{},{:.3},{}",
        pu,
        metrics.lengths.len(),
        metrics.yield_bases,
        calculations::get_n(&metrics.lengths, metrics.yield_bases, 0.50),
        calculations::mean_accuracy(metrics.qual_sum, metrics.qual_reads),
        metrics.modifications
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::{Header, HeaderRecord};

    #[test]
    fn verify_app() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }

    #[test]
    fn extract_pu_joins_unique_pus_and_falls_back_cleanly() {
        let empty = Header::new();
        assert_eq!(extract_pu(&empty, "dir/sample.bam"), "sample.bam");
        assert_eq!(extract_pu(&empty, "-"), "stdin");

        let mut no_pu = Header::new();
        let mut rg1 = HeaderRecord::new(b"RG");
        rg1.push_tag(b"ID", "rg1");
        let mut rg2 = HeaderRecord::new(b"RG");
        rg2.push_tag(b"ID", "rg2");
        let mut rg3 = HeaderRecord::new(b"RG");
        rg3.push_tag(b"ID", "rg1");
        no_pu.push_record(&rg1).push_record(&rg2).push_record(&rg3);
        assert_eq!(extract_pu(&no_pu, "dir/sample.bam"), "rg1;rg2");

        let mut multi_pu = Header::new();
        let mut r1 = HeaderRecord::new(b"RG");
        r1.push_tag(b"ID", "rg1");
        r1.push_tag(b"PU", "fc_1");
        let mut r2 = HeaderRecord::new(b"RG");
        r2.push_tag(b"ID", "rg2");
        r2.push_tag(b"PU", "fc_2");
        let mut r3 = HeaderRecord::new(b"RG");
        r3.push_tag(b"ID", "rg3");
        r3.push_tag(b"PU", "fc_1");
        multi_pu.push_record(&r1).push_record(&r2).push_record(&r3);
        assert_eq!(extract_pu(&multi_pu, "dir/sample.bam"), "fc_1;fc_2");
    }
}
