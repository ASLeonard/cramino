use rayon::prelude::*;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use rust_htslib::{bam, htslib};
use std::collections::BTreeSet;
use std::sync::LazyLock;

static PHRED_TO_ERROR: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0_f64; 256];
    for (q, slot) in table.iter_mut().enumerate() {
        *slot = 10.0_f64.powf(-(q as f64) / 10.0);
    }
    table
});

const MOD_5MC: u32 = 1 << 0;
const MOD_5HMC: u32 = 1 << 1;
const MOD_5FC: u32 = 1 << 2;
const MOD_5CAC: u32 = 1 << 3;
const MOD_4MC: u32 = 1 << 4;
const MOD_6MA: u32 = 1 << 5;
const MOD_8OXOG: u32 = 1 << 6;
const MOD_5HMU: u32 = 1 << 7;
const MOD_5FU: u32 = 1 << 8;
const MOD_5CAU: u32 = 1 << 9;
const MOD_PSEU: u32 = 1 << 10;
const MOD_INOSINE: u32 = 1 << 11;
const MOD_2OME_A: u32 = 1 << 12;
const MOD_2OME_C: u32 = 1 << 13;
const MOD_2OME_G: u32 = 1 << 14;
const MOD_2OME_U: u32 = 1 << 15;
const MOD_XAO: u32 = 1 << 16;
const MOD_C: u32 = 1 << 17;
const MOD_A: u32 = 1 << 18;
const MOD_G: u32 = 1 << 19;
const MOD_T: u32 = 1 << 20;
const MOD_U: u32 = 1 << 21;
const MOD_N: u32 = 1 << 22;

const KNOWN_MODS: [(u32, &str); 23] = [
    (MOD_5MC, "5mC"),
    (MOD_5HMC, "5hmC"),
    (MOD_5FC, "5fC"),
    (MOD_5CAC, "5caC"),
    (MOD_4MC, "4mC"),
    (MOD_6MA, "6mA"),
    (MOD_8OXOG, "8oxoG"),
    (MOD_5HMU, "5hmU"),
    (MOD_5FU, "5fU"),
    (MOD_5CAU, "5caU"),
    (MOD_PSEU, "pseU"),
    (MOD_INOSINE, "inosine"),
    (MOD_2OME_A, "2'-O-Me-A"),
    (MOD_2OME_C, "2'-O-Me-C"),
    (MOD_2OME_G, "2'-O-Me-G"),
    (MOD_2OME_U, "2'-O-Me-U"),
    (MOD_XAO, "Xao"),
    (MOD_C, "modC"),
    (MOD_A, "modA"),
    (MOD_G, "modG"),
    (MOD_T, "modT"),
    (MOD_U, "modU"),
    (MOD_N, "modN"),
];

pub struct Data {
    pub lengths: Vec<u32>,
    pub yield_bases: u64,
    pub qual_sum: f64,
    pub qual_reads: usize,
    pub modifications: String,
}

#[derive(Default)]
struct ModSummary {
    known_mask: u32,
    unknown: BTreeSet<String>,
}

impl ModSummary {
    fn merge(&mut self, known_mask: u32, unknown: Vec<String>) {
        self.known_mask |= known_mask;
        if !unknown.is_empty() {
            self.unknown.extend(unknown);
        }
    }

    fn format(&self) -> String {
        if self.known_mask == 0 && self.unknown.is_empty() {
            return "None".to_string();
        }
        let mut names = Vec::new();
        for &(bit, name) in &KNOWN_MODS {
            if (self.known_mask & bit) != 0 {
                names.push(name);
            }
        }
        for name in &self.unknown {
            names.push(name.as_str());
        }
        names.join(";")
    }
}

/// Extract read lengths, total yield, quality error-rate sum, and base modification types.
pub fn extract(args: &crate::Cli) -> (Data, rust_htslib::bam::Header) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    let mut bam = if args.input == "-" {
        bam::Reader::from_stdin().expect(
            "\n\nError reading alignments from stdin.\nDid you include the file header with -h?\n\n\n\n",
        )
    } else {
        bam::Reader::from_path(&args.input)
            .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n")
    };

    if args.input.ends_with(".cram") && args.reference.is_some() {
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
            hts_sys::sam_fields_SAM_FLAG
                | hts_sys::sam_fields_SAM_SEQ
                | hts_sys::sam_fields_SAM_QUAL
                | hts_sys::sam_fields_SAM_AUX,
        )
        .expect("Failed setting cram options");
    }

    let header = rust_htslib::bam::Header::from_template(bam.header());
    if args.threads > 0 {
        bam.set_threads(args.threads)
            .expect("Failure setting decompression threads");
    }

    let filter_mask = (htslib::BAM_FSECONDARY | htslib::BAM_FSUPPLEMENTARY) as u16;
    let mut lengths = vec![];
    let mut yield_bases: u64 = 0;
    let mut qual_sum: f64 = 0.0;
    let mut qual_reads: usize = 0;
    let mut mod_summary = ModSummary::default();
    let mut record = bam::Record::new();

    while let Some(res) = bam.read(&mut record) {
        res.expect("Failure parsing Bam file");
        if (record.flags() & filter_mask) != 0 {
            continue;
        }
        let len = record.seq_len() as u32;
        lengths.push(len);
        yield_bases += u64::from(len);

        if let Some(err) = get_read_error_rate(&record) {
            qual_sum += err;
            qual_reads += 1;
        }
        let (mask, unknown) = get_mod_tags(&record);
        mod_summary.merge(mask, unknown);
    }

    lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    (
        Data {
            lengths,
            yield_bases,
            qual_sum,
            qual_reads,
            modifications: mod_summary.format(),
        },
        header,
    )
}

fn get_read_error_rate(record: &bam::Record) -> Option<f64> {
    if let Ok(aux) = record.aux(b"rq") {
        let rq = match aux {
            Aux::Float(v) => Some(f64::from(v)),
            Aux::Double(v) => Some(v),
            _ => None,
        };
        if let Some(v) = rq.filter(|v| v.is_finite() && (0.0..=1.0).contains(v)) {
            return Some(1.0 - v);
        }
    }

    if let Ok(aux) = record.aux(b"qs") {
        let qs = match aux {
            Aux::Float(v) => Some(f64::from(v)),
            Aux::Double(v) => Some(v),
            Aux::I8(v) => Some(f64::from(v)),
            Aux::U8(v) => Some(f64::from(v)),
            Aux::I16(v) => Some(f64::from(v)),
            Aux::U16(v) => Some(f64::from(v)),
            Aux::I32(v) => Some(f64::from(v)),
            Aux::U32(v) => Some(f64::from(v)),
            _ => None,
        };
        if let Some(v) = qs.filter(|v| v.is_finite() && *v >= 0.0) {
            return Some(10.0_f64.powf(-v / 10.0));
        }
    }

    mean_qual_error(record.qual())
}

#[inline]
fn mean_qual_error(qual: &[u8]) -> Option<f64> {
    if qual.is_empty() || qual[0] == 0xFF {
        return None;
    }
    let table = &*PHRED_TO_ERROR;
    let mut sum = 0.0_f64;
    for &q in qual {
        sum += table[q as usize];
    }
    Some(sum / qual.len() as f64)
}

fn get_mod_tags(record: &bam::Record) -> (u32, Vec<String>) {
    let aux = record.aux(b"MM").or_else(|_| record.aux(b"Mm"));
    match aux {
        Ok(Aux::String(mm)) => parse_mm_bytes(mm.as_bytes()),
        _ => (0, Vec::new()),
    }
}

fn parse_mm_bytes(mm: &[u8]) -> (u32, Vec<String>) {
    let mut known_mask = 0u32;
    let mut unknown = Vec::new();
    let mut pos = 0;

    while pos < mm.len() {
        let series_end = mm[pos..]
            .iter()
            .position(|&b| b == b';')
            .map_or(mm.len(), |rel| pos + rel);
        let series = &mm[pos..series_end];
        pos = series_end + 1;

        if series.is_empty() {
            continue;
        }

        // Only the header prefix before the first ',' contains base, strand, and mod code(s).
        let header_end = series
            .iter()
            .position(|&b| b == b',')
            .unwrap_or(series.len());
        let mut header = &series[..header_end];
        if let Some(&last) = header.last() {
            if last == b'.' || last == b'?' {
                header = &header[..header.len() - 1];
            }
        }
        if header.len() < 3
            || !matches!(header[0], b'A' | b'C' | b'G' | b'T' | b'U' | b'N')
            || !matches!(header[1], b'+' | b'-')
        {
            continue;
        }

        let base = header[0] as char;
        let strand = header[1] as char;
        let mod_code = &header[2..];

        if mod_code.iter().all(|b| b.is_ascii_digit()) {
            let code_str = std::str::from_utf8(mod_code).unwrap();
            if let Some(bit) = chebi_mod_mask(code_str) {
                known_mask |= bit;
            } else {
                unknown.push(format!("{base}{strand}{code_str}"));
            }
        } else if mod_code.iter().all(|b| b.is_ascii_alphabetic()) {
            for &code_byte in mod_code {
                if let Some(bit) = single_char_mod_mask(code_byte) {
                    known_mask |= bit;
                } else {
                    unknown.push(format!("{base}{strand}{}", code_byte as char));
                }
            }
        }
    }

    (known_mask, unknown)
}

fn single_char_mod_mask(code: u8) -> Option<u32> {
    match code {
        b'm' => Some(MOD_5MC),
        b'h' => Some(MOD_5HMC),
        b'f' => Some(MOD_5FC),
        b'c' => Some(MOD_5CAC),
        b'a' => Some(MOD_6MA),
        b'o' => Some(MOD_8OXOG),
        b'g' => Some(MOD_5HMU),
        b'e' => Some(MOD_5FU),
        b'b' => Some(MOD_5CAU),
        b'n' => Some(MOD_XAO),
        b'C' => Some(MOD_C),
        b'A' => Some(MOD_A),
        b'G' => Some(MOD_G),
        b'T' => Some(MOD_T),
        b'U' => Some(MOD_U),
        b'N' => Some(MOD_N),
        _ => None,
    }
}

fn chebi_mod_mask(code: &str) -> Option<u32> {
    match code {
        "27551" => Some(MOD_5MC),
        "76792" => Some(MOD_5HMC),
        "76794" => Some(MOD_5FC),
        "76793" => Some(MOD_5CAC),
        "21839" => Some(MOD_4MC),
        "28871" => Some(MOD_6MA),
        "44605" => Some(MOD_8OXOG),
        "76795" => Some(MOD_5HMU),
        "80961" => Some(MOD_5FU),
        "76796" => Some(MOD_5CAU),
        "17802" => Some(MOD_PSEU),
        "17596" => Some(MOD_INOSINE),
        "69426" => Some(MOD_2OME_A),
        "19228" => Some(MOD_2OME_C),
        "19229" => Some(MOD_2OME_G),
        "19227" => Some(MOD_2OME_U),
        "18107" => Some(MOD_XAO),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::{Header, HeaderRecord};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn format_mm(mm: &str) -> String {
        let mut summary = ModSummary::default();
        let (mask, unknown) = parse_mm_bytes(mm.as_bytes());
        summary.merge(mask, unknown);
        summary.format()
    }

    fn make_record(seq: &[u8], rq: Option<f32>, mm_tag: Option<(&[u8; 2], &str)>) -> bam::Record {
        let mut record = bam::Record::new();
        let qual = vec![30u8; seq.len()];
        record.set(b"read1", None, seq, &qual);
        if let Some(v) = rq {
            record.push_aux(b"rq", Aux::Float(v)).unwrap();
        }
        if let Some((tag, val)) = mm_tag {
            record.push_aux(tag, Aux::String(val)).unwrap();
        }
        record
    }

    #[test]
    fn parse_standard_and_chebi_mm_tags() {
        assert_eq!(format_mm("C+m?,0,2,5;"), "5mC");
        assert_eq!(format_mm("C+h?,0,1;C+m?,2,3;"), "5mC;5hmC");
        assert_eq!(format_mm("C+mh?,0,1;A+a.,0;"), "5mC;5hmC;6mA");
        assert_eq!(format_mm("C+21839?,0,1;C+m?,0;"), "5mC;4mC");
        assert_eq!(format_mm("A+17596.;T+17802?;"), "pseU;inosine");
        assert_eq!(format_mm("C+m?;"), "5mC");
        assert_eq!(format_mm(""), "None");
    }

    #[test]
    fn parse_ambiguous_and_unknown_mm_tags() {
        assert_eq!(format_mm("C+C?,0;A+A.,1;"), "modC;modA");
        assert_eq!(format_mm("C+99999?,0;C+m?,1;"), "5mC;C+99999");
        assert_eq!(format_mm("G-99999?,0;"), "G-99999");
    }

    #[test]
    fn read_error_rate_and_mod_tags_handle_rq_qs_qual_and_legacy_mm() {
        let mut rec_mm = make_record(b"ACGTACGT", Some(0.99), Some((b"MM", "C+mh?,0;")));
        rec_mm
            .push_aux(b"ML", Aux::ArrayU8((&[250u8, 12u8][..]).into()))
            .unwrap();
        assert!((get_read_error_rate(&rec_mm).unwrap() - 0.01).abs() < 1e-6);
        assert_eq!(get_mod_tags(&rec_mm), (MOD_5MC | MOD_5HMC, Vec::new()));

        let mut rec_qs_float = make_record(b"ACGTACGT", None, Some((b"Mm", "A+a.,0;")));
        rec_qs_float.push_aux(b"qs", Aux::Float(20.0)).unwrap();
        assert!((get_read_error_rate(&rec_qs_float).unwrap() - 0.01).abs() < 1e-6);
        assert_eq!(get_mod_tags(&rec_qs_float), (MOD_6MA, Vec::new()));

        let mut rec_qs_then_rq = make_record(b"ACGTACGT", None, None);
        rec_qs_then_rq.push_aux(b"qs", Aux::I32(10)).unwrap();
        rec_qs_then_rq.push_aux(b"rq", Aux::Float(0.99)).unwrap();
        assert!((get_read_error_rate(&rec_qs_then_rq).unwrap() - 0.01).abs() < 1e-6);

        let rec_none = make_record(b"ACGTACGT", None, None);
        assert_eq!(get_mod_tags(&rec_none), (0, Vec::new()));
        assert!((get_read_error_rate(&rec_none).unwrap() - 0.001).abs() < 1e-9);
        assert_eq!(mean_qual_error(&[0xFF, 0xFF, 0xFF]), None);
    }

    #[test]
    fn extract_handles_synthetic_and_repo_bams() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let bam_path = std::env::temp_dir().join(format!("cramino_test_{nanos}.bam"));

        let mut header = Header::new();
        let mut rg = HeaderRecord::new(b"RG");
        rg.push_tag(b"ID", "rg1");
        rg.push_tag(b"PU", "flowcell_pu");
        header.push_record(&rg);

        {
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut r1 = make_record(b"ACGTACGTAC", Some(0.99), Some((b"MM", "C+mh?,0;")));
            r1.push_aux(b"ML", Aux::ArrayU8((&[240u8, 10u8][..]).into()))
                .unwrap();
            let mut r2 = make_record(b"AAACCCGGGTTT", None, Some((b"Mm", "A+a.,0;C+21839?;")));
            r2.push_aux(b"qs", Aux::Float(30.0)).unwrap();

            let mut r3 = bam::Record::new();
            r3.set(b"illumina_read", None, b"ACGT", &[20, 20, 20, 20]);

            let mut r4 = bam::Record::new();
            r4.set(b"unscored", None, b"ACGTAC", &[0xFF; 6]);
            r4.push_aux(b"rq", Aux::Float(-1.0)).unwrap();

            let mut r_sec =
                make_record(b"TTTTTTTTTTTTTTTTTTTT", Some(0.5), Some((b"MM", "G+o?,0;")));
            r_sec.set_flags(0x100);
            let mut r_sup =
                make_record(b"CCCCCCCCCCCCCCCCCCCC", Some(0.5), Some((b"MM", "T+g?,0;")));
            r_sup.set_flags(0x800);

            writer.write(&r1).unwrap();
            writer.write(&r2).unwrap();
            writer.write(&r3).unwrap();
            writer.write(&r4).unwrap();
            writer.write(&r_sec).unwrap();
            writer.write(&r_sup).unwrap();
        }

        let cli = crate::Cli {
            input: bam_path.to_str().unwrap().to_string(),
            threads: 2,
            reference: None,
        };
        let (data, parsed_header) = extract(&cli);
        std::fs::remove_file(&bam_path).ok();

        assert_eq!(data.lengths, vec![12, 10, 6, 4]);
        assert_eq!(data.yield_bases, 32);
        assert_eq!(data.qual_reads, 3);
        assert!((data.qual_sum - (0.01 + 0.001 + 0.01)).abs() < 1e-6);
        assert_eq!(data.modifications, "5mC;5hmC;4mC;6mA");
        assert_eq!(parsed_header.to_hashmap()["RG"][0]["PU"], "flowcell_pu");

        let phased_cli = crate::Cli {
            input: "test-data/small-test-phased.bam".to_string(),
            threads: 4,
            reference: None,
        };
        let ubam_cli = crate::Cli {
            input: "test-data/small-test-ubam.bam".to_string(),
            threads: 4,
            reference: None,
        };
        let (phased_data, _) = extract(&phased_cli);
        let (ubam_data, _) = extract(&ubam_cli);
        assert_eq!(phased_data.lengths, ubam_data.lengths);
        assert_eq!(phased_data.yield_bases, ubam_data.yield_bases);
        assert_eq!(phased_data.qual_reads, ubam_data.qual_reads);
        assert!((phased_data.qual_sum - ubam_data.qual_sum).abs() < 1e-9);
        assert_eq!(phased_data.modifications, "None");
    }
}
