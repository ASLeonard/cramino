use rust_htslib::bam::record::{Aux};
use rust_htslib::{bam, bam::Read};

pub struct Data {
    pub lengths: Option<Vec<u128>>,
    pub all_counts: usize,
    pub qualities: Option<Vec<f64>>,
}

pub fn extract(args: &crate::Cli) -> (Data, rust_htslib::bam::Header) {
    let mut lengths = vec![];
    let mut qualities = vec![];
    let mut bam = if args.input == "-" {
        bam::Reader::from_stdin().expect("\n\nError reading alignments from stdin.\nDid you include the file header with -h?\n\n\n\n")
    } else {
        bam::Reader::from_path(&args.input)
            .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n")
    };
    if args.input.ends_with(".cram") & args.reference.is_some() {
        // bam.set_cram_option(htslib::CFR_REQUIRED_FIELDS, htslib::sam_fields_SAM_AUX as i32)
        //     .expect("Failed setting cram options");
        bam.set_reference(
            args.reference
                .as_ref()
                .expect("Failed setting reference for CRAM file"),
        )
        .expect("Failed setting reference for CRAM file");
    }
    if args.input.ends_with(".cram") {
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_CIGAR
                | hts_sys::sam_fields_SAM_SEQ,
        )
        .expect("Failed setting cram options");
    }
    let header = bam.header().clone();
    let header = rust_htslib::bam::Header::from_template(&header);
    bam.set_threads(args.threads)
        .expect("Failure setting decompression threads");

    let mut all_counts = 0;
    for read in bam
        .rc_records()
        .map(|r| r.expect("Failure parsing Bam file"))
        .inspect(|_| all_counts += 1)
    {
        lengths.push(read.seq_len() as u128);
        qualities.push(get_rq_tag(&read).into());
    }
    // sort vectors in descending order (required for N50/N75)
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    (
        Data {
            lengths: Some(lengths),
            all_counts,
            qualities: Some(qualities),
            },
        header,
    )
}



fn get_rq_tag(record: &bam::Record) -> f32 {
    match record.aux(b"rq") {
        Ok(value) => match value {
            Aux::Float(v) =>
                v,
		_ => todo!()
        },
        Err(_e) => -1.0,
    }
}
