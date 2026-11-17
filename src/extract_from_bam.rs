use rayon::prelude::*;
use rust_htslib::bam::Read;
use rust_htslib::{bam, htslib};
use std::collections::{BTreeSet, VecDeque};
use std::fs::File;
#[cfg(unix)]
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::sync::LazyLock;

const BGZF_MAX_BLOCK_SIZE: usize = 65_536;
const COMP_READ_BUF_SIZE: usize = 1024 * 1024;
const STREAM_COMPACT_THRESHOLD: usize = 256 * 1024;
const MAX_BAM_RECORD_LEN: usize = 256 * 1024 * 1024;
const MAX_BAM_SEQ_LEN: i32 = 100_000_000;
const BAM_SECONDARY_OR_SUPPLEMENTARY: u16 = 0x100 | 0x800;

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

mod libdeflate {
    use std::ffi::c_void;

    const LIBDEFLATE_SUCCESS: u32 = 0;

    #[repr(C)]
    struct LibdeflateDecompressor {
        _private: [u8; 0],
    }

    extern "C" {
        fn libdeflate_alloc_decompressor() -> *mut LibdeflateDecompressor;
        fn libdeflate_free_decompressor(decompressor: *mut LibdeflateDecompressor);
        fn libdeflate_deflate_decompress(
            decompressor: *mut LibdeflateDecompressor,
            in_buf: *const c_void,
            in_nbytes: usize,
            out_buf: *mut c_void,
            out_nbytes_avail: usize,
            actual_out_nbytes_ret: *mut usize,
        ) -> u32;
    }

    pub struct Decompressor {
        ptr: *mut LibdeflateDecompressor,
    }

    impl Decompressor {
        pub fn new() -> Self {
            let ptr = unsafe { libdeflate_alloc_decompressor() };
            assert!(!ptr.is_null(), "Failed to allocate libdeflate decompressor");
            Self { ptr }
        }

        pub fn decompress_deflate(
            &mut self,
            input: &[u8],
            output: &mut [u8],
            expected_isize: usize,
        ) -> Option<usize> {
            let mut actual_out = 0usize;
            let rc = unsafe {
                libdeflate_deflate_decompress(
                    self.ptr,
                    input.as_ptr().cast(),
                    input.len(),
                    output.as_mut_ptr().cast(),
                    output.len(),
                    &mut actual_out,
                )
            };
            if rc == LIBDEFLATE_SUCCESS && actual_out == expected_isize {
                Some(actual_out)
            } else {
                None
            }
        }
    }

    impl Drop for Decompressor {
        fn drop(&mut self) {
            unsafe { libdeflate_free_decompressor(self.ptr) };
        }
    }
}

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

    fn merge_summary(&mut self, other: ModSummary) {
        self.known_mask |= other.known_mask;
        if !other.unknown.is_empty() {
            self.unknown.extend(other.unknown);
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

#[derive(Default)]
struct ChunkOutput {
    lengths: Vec<u32>,
    yield_bases: u64,
    qual_sum: f64,
    qual_reads: usize,
    mods: ModSummary,
}

#[cfg(unix)]
struct CompressedReader<'a> {
    file: &'a File,
    file_size: u64,
    comp_buf: Vec<u8>,
    comp_start_off: u64,
    comp_pos: usize,
    comp_len: usize,
    next_block_off: u64,
    decompressor: libdeflate::Decompressor,
    eof: bool,
}

#[cfg(unix)]
impl<'a> CompressedReader<'a> {
    fn new(file: &'a File, file_size: u64, start_block_off: u64) -> Self {
        Self {
            file,
            file_size,
            comp_buf: vec![0u8; COMP_READ_BUF_SIZE],
            comp_start_off: start_block_off,
            comp_pos: 0,
            comp_len: 0,
            next_block_off: start_block_off,
            decompressor: libdeflate::Decompressor::new(),
            eof: false,
        }
    }

    fn ensure_comp_bytes(&mut self, needed: usize) -> bool {
        if self.comp_len - self.comp_pos >= needed {
            return true;
        }
        if needed > self.comp_buf.len() {
            self.comp_buf.resize(needed, 0);
        }

        let rem = self.comp_len - self.comp_pos;
        if rem > 0 && self.comp_pos > 0 {
            self.comp_buf.copy_within(self.comp_pos..self.comp_len, 0);
        }
        self.comp_start_off = self.next_block_off;
        self.comp_pos = 0;
        self.comp_len = rem;

        while self.comp_len < needed {
            let read_off = self.comp_start_off + self.comp_len as u64;
            if read_off >= self.file_size {
                break;
            }
            let max_to_read =
                (self.file_size - read_off).min((self.comp_buf.len() - self.comp_len) as u64) as usize;
            let n = self
                .file
                .read_at(&mut self.comp_buf[self.comp_len..self.comp_len + max_to_read], read_off)
                .expect("Failure reading BAM file");
            if n == 0 {
                break;
            }
            self.comp_len += n;
        }

        self.comp_len - self.comp_pos >= needed
    }

    /// Reads the next BGZF block header and footer metadata `(bsize, header_len, isize)` without decompressing.
    #[inline]
    fn peek_next_block_meta(&mut self) -> Option<(usize, usize, usize)> {
        if self.eof || self.next_block_off >= self.file_size {
            self.eof = true;
            return None;
        }
        if !self.ensure_comp_bytes(18) {
            panic!("Truncated BGZF block header at offset {}", self.next_block_off);
        }
        let (bsize, header_len) =
            parse_bgzf_header(&self.comp_buf[self.comp_pos..self.comp_pos + 18])
                .unwrap_or_else(|| panic!("Invalid BGZF block header at offset {}", self.next_block_off));
        if !self.ensure_comp_bytes(bsize) {
            panic!("Truncated BGZF block at offset {}", self.next_block_off);
        }
        let block = &self.comp_buf[self.comp_pos..self.comp_pos + bsize];
        let isize = u32::from_le_bytes([
            block[bsize - 4],
            block[bsize - 3],
            block[bsize - 2],
            block[bsize - 1],
        ]) as usize;
        assert!(
            isize <= BGZF_MAX_BLOCK_SIZE,
            "Invalid BGZF uncompressed size {isize} at offset {}",
            self.next_block_off
        );
        Some((bsize, header_len, isize))
    }

    #[inline]
    fn consume_compressed_block(&mut self, bsize: usize) -> u64 {
        let block_file_off = self.next_block_off;
        self.comp_pos += bsize;
        self.next_block_off += bsize as u64;
        block_file_off
    }

    #[inline]
    fn decompress_current_block(
        &mut self,
        bsize: usize,
        header_len: usize,
        isize: usize,
        out: &mut [u8],
    ) -> u64 {
        let block_file_off = self.next_block_off;
        if isize > 0 {
            let deflate_slice = &self.comp_buf[self.comp_pos + header_len..self.comp_pos + bsize - 8];
            self.decompressor
                .decompress_deflate(deflate_slice, out, isize)
                .unwrap_or_else(|| panic!("Corrupted BGZF block at offset {block_file_off}"));
        }
        self.comp_pos += bsize;
        self.next_block_off += bsize as u64;
        block_file_off
    }
}

struct StreamBuf {
    data: Vec<u8>,
    cursor: usize,
    /// `(start_in_data, compressed_file_offset)` for non-empty BGZF blocks in `data`.
    blocks: VecDeque<(usize, u64)>,
}

impl StreamBuf {
    fn new() -> Self {
        Self {
            data: Vec::with_capacity(512 * 1024),
            cursor: 0,
            blocks: VecDeque::with_capacity(16),
        }
    }

    #[inline]
    fn block_file_offset_at(&self, pos: usize) -> Option<u64> {
        let idx = self.blocks.partition_point(|&(start, _)| start <= pos);
        if idx == 0 {
            None
        } else {
            Some(self.blocks[idx - 1].1)
        }
    }

    fn compact(&mut self) {
        if self.cursor == 0 {
            return;
        }
        if self.cursor >= self.data.len() {
            self.data.clear();
            self.blocks.clear();
            self.cursor = 0;
            return;
        }

        let keep_idx = self.blocks.partition_point(|&(start, _)| start <= self.cursor);
        if keep_idx > 1 {
            self.blocks.drain(..keep_idx - 1);
        }
        let remaining = self.data.len() - self.cursor;
        self.data.copy_within(self.cursor.., 0);
        self.data.truncate(remaining);

        for (i, (start, _)) in self.blocks.iter_mut().enumerate() {
            if i == 0 && *start <= self.cursor {
                *start = 0;
            } else {
                *start -= self.cursor;
            }
        }
        self.cursor = 0;
    }
}

#[cfg(unix)]
struct SyncStream<'a> {
    comp: CompressedReader<'a>,
    stream: StreamBuf,
}

#[cfg(unix)]
impl<'a> SyncStream<'a> {
    fn new(file: &'a File, file_size: u64, start_block_off: u64) -> Self {
        Self {
            comp: CompressedReader::new(file, file_size, start_block_off),
            stream: StreamBuf::new(),
        }
    }

    fn next_bgzf_block(&mut self) -> bool {
        while let Some((bsize, header_len, isize)) = self.comp.peek_next_block_meta() {
            if isize == 0 {
                self.comp.consume_compressed_block(bsize);
                continue;
            }
            let old_len = self.stream.data.len();
            self.stream.data.reserve(isize);
            // SAFETY: `reserve(isize)` ensures capacity >= `old_len + isize` and `decompress_deflate`
            // returns `Some(isize)` only after writing all `isize` bytes.
            let out_slice = unsafe {
                std::slice::from_raw_parts_mut(self.stream.data.as_mut_ptr().add(old_len), isize)
            };
            let block_file_off =
                self.comp
                    .decompress_current_block(bsize, header_len, isize, out_slice);
            self.stream.blocks.push_back((old_len, block_file_off));
            unsafe {
                self.stream.data.set_len(old_len + isize);
            }
            return true;
        }
        false
    }

    #[inline]
    fn ensure_uncomp_bytes(&mut self, required_end: usize) -> bool {
        while self.stream.data.len() < required_end {
            if !self.next_bgzf_block() {
                return false;
            }
        }
        true
    }

    fn into_block_reader(self, p: usize) -> BlockReader<'a> {
        let b_idx = self
            .stream
            .blocks
            .partition_point(|&(start, _)| start <= p)
            .saturating_sub(1);
        let (b_start, b_off) = self.stream.blocks[b_idx];
        let b_end = self
            .stream
            .blocks
            .get(b_idx + 1)
            .map_or(self.stream.data.len(), |b| b.0);
        let b_len = b_end - b_start;

        let mut block_buf = vec![0u8; BGZF_MAX_BLOCK_SIZE];
        block_buf[..b_len].copy_from_slice(&self.stream.data[b_start..b_end]);

        let mut queued_blocks = VecDeque::new();
        for j in (b_idx + 1)..self.stream.blocks.len() {
            let (j_start, j_off) = self.stream.blocks[j];
            let j_end = self
                .stream
                .blocks
                .get(j + 1)
                .map_or(self.stream.data.len(), |b| b.0);
            queued_blocks.push_back((j_off, self.stream.data[j_start..j_end].to_vec()));
        }

        BlockReader {
            comp: self.comp,
            block_buf,
            block_pos: p - b_start,
            block_len: b_len,
            curr_block_file_off: b_off,
            queued_blocks,
        }
    }
}

#[cfg(unix)]
struct BlockReader<'a> {
    comp: CompressedReader<'a>,
    block_buf: Vec<u8>,
    block_pos: usize,
    block_len: usize,
    curr_block_file_off: u64,
    queued_blocks: VecDeque<(u64, Vec<u8>)>,
}

#[cfg(unix)]
impl<'a> BlockReader<'a> {
    fn new(file: &'a File, file_size: u64, start_block_off: u64) -> Self {
        Self {
            comp: CompressedReader::new(file, file_size, start_block_off),
            block_buf: vec![0u8; BGZF_MAX_BLOCK_SIZE],
            block_pos: 0,
            block_len: 0,
            curr_block_file_off: start_block_off,
            queued_blocks: VecDeque::new(),
        }
    }

    #[inline]
    fn avail(&self) -> usize {
        self.block_len - self.block_pos
    }

    #[inline]
    fn ensure_nonempty_block(&mut self) -> bool {
        if self.block_pos < self.block_len {
            return true;
        }
        self.load_next_block()
    }

    fn load_next_block(&mut self) -> bool {
        if let Some((off, qb)) = self.queued_blocks.pop_front() {
            let len = qb.len();
            self.block_buf[..len].copy_from_slice(&qb);
            self.block_pos = 0;
            self.block_len = len;
            self.curr_block_file_off = off;
            return true;
        }
        while let Some((bsize, header_len, isize)) = self.comp.peek_next_block_meta() {
            if isize == 0 {
                self.comp.consume_compressed_block(bsize);
                continue;
            }
            let off = self.comp.decompress_current_block(
                bsize,
                header_len,
                isize,
                &mut self.block_buf[..isize],
            );
            self.block_pos = 0;
            self.block_len = isize;
            self.curr_block_file_off = off;
            return true;
        }
        self.block_pos = 0;
        self.block_len = 0;
        false
    }

    #[inline]
    fn read_exact(&mut self, mut dst: &mut [u8]) -> bool {
        if self.block_pos + dst.len() <= self.block_len {
            dst.copy_from_slice(&self.block_buf[self.block_pos..self.block_pos + dst.len()]);
            self.block_pos += dst.len();
            return true;
        }
        while !dst.is_empty() {
            if !self.ensure_nonempty_block() {
                return false;
            }
            let take = dst.len().min(self.block_len - self.block_pos);
            dst[..take].copy_from_slice(&self.block_buf[self.block_pos..self.block_pos + take]);
            self.block_pos += take;
            dst = &mut dst[take..];
        }
        true
    }

    /// Advances by `n` uncompressed bytes, skipping whole BGZF blocks without decompression.
    #[inline]
    fn skip_uncomp_bytes(&mut self, mut n: usize) -> bool {
        let avail = self.block_len - self.block_pos;
        if n <= avail {
            self.block_pos += n;
            return true;
        }
        n -= avail;
        self.block_pos = self.block_len;

        while let Some((_, qb)) = self.queued_blocks.front() {
            if n >= qb.len() {
                let (off, qb) = self.queued_blocks.pop_front().unwrap();
                n -= qb.len();
                self.curr_block_file_off = off;
                if n == 0 {
                    return true;
                }
            } else {
                let (off, qb) = self.queued_blocks.pop_front().unwrap();
                let len = qb.len();
                self.block_buf[..len].copy_from_slice(&qb);
                self.block_len = len;
                self.block_pos = n;
                self.curr_block_file_off = off;
                return true;
            }
        }

        while n > 0 {
            let Some((bsize, header_len, isize)) = self.comp.peek_next_block_meta() else {
                return false;
            };
            if isize == 0 {
                self.comp.consume_compressed_block(bsize);
                continue;
            }
            if n >= isize {
                self.curr_block_file_off = self.comp.consume_compressed_block(bsize);
                n -= isize;
                if n == 0 {
                    self.block_pos = 0;
                    self.block_len = 0;
                    return true;
                }
            } else {
                self.curr_block_file_off = self.comp.decompress_current_block(
                    bsize,
                    header_len,
                    isize,
                    &mut self.block_buf[..isize],
                );
                self.block_len = isize;
                self.block_pos = n;
                return true;
            }
        }

        true
    }

    fn skip_until_nul(&mut self, max_len: usize) -> usize {
        let mut consumed = 0;
        while consumed < max_len {
            if !self.ensure_nonempty_block() {
                break;
            }
            let avail = (self.block_len - self.block_pos).min(max_len - consumed);
            let slice = &self.block_buf[self.block_pos..self.block_pos + avail];
            if let Some(pos) = slice.iter().position(|&b| b == 0) {
                let take = pos + 1;
                self.block_pos += take;
                return consumed + take;
            }
            self.block_pos += avail;
            consumed += avail;
        }
        consumed
    }

    fn read_z_into(&mut self, dst: &mut Vec<u8>, max_len: usize) -> usize {
        dst.clear();
        let mut consumed = 0;
        while consumed < max_len {
            if !self.ensure_nonempty_block() {
                break;
            }
            let avail = (self.block_len - self.block_pos).min(max_len - consumed);
            let slice = &self.block_buf[self.block_pos..self.block_pos + avail];
            if let Some(pos) = slice.iter().position(|&b| b == 0) {
                dst.extend_from_slice(&slice[..pos]);
                let take = pos + 1;
                self.block_pos += take;
                return consumed + take;
            }
            dst.extend_from_slice(slice);
            self.block_pos += avail;
            consumed += avail;
        }
        consumed
    }
}

/// Extract read lengths, total yield, quality error-rate sum, and base modification types.
pub fn extract(args: &crate::Cli) -> (Data, rust_htslib::bam::Header) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    let bam = if args.input == "-" {
        bam::Reader::from_stdin().expect(
            "\n\nError reading alignments from stdin.\nDid you include the file header with -h?\n\n\n\n",
        )
    } else {
        bam::Reader::from_path(&args.input)
            .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n")
    };

    let header = rust_htslib::bam::Header::from_template(bam.header());
    let n_targets = bam.header().target_count() as i32;

    #[cfg(unix)]
    if can_read_parallel_bam(args, &bam) {
        let voffset = bam.tell();
        drop(bam);
        let data = extract_parallel_bam(&args.input, args.threads, n_targets, voffset);
        return (data, header);
    }

    let data = extract_streaming(args, bam);
    (data, header)
}

#[cfg(unix)]
fn can_read_parallel_bam(args: &crate::Cli, bam: &bam::Reader) -> bool {
    if args.input == "-" || args.input.ends_with(".cram") {
        return false;
    }
    let fmt = unsafe { &(*bam.htsfile()).format };
    if fmt.format != htslib::htsExactFormat_bam || fmt.compression != htslib::htsCompression_bgzf {
        return false;
    }
    std::fs::metadata(&args.input)
        .map(|m| m.is_file())
        .unwrap_or(false)
}

fn extract_streaming(args: &crate::Cli, mut bam: bam::Reader) -> Data {
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
    if args.threads > 0 {
        bam.set_threads(args.threads)
            .expect("Failure setting decompression threads");
    }

    let mut lengths = vec![];
    let mut yield_bases: u64 = 0;
    let mut qual_sum: f64 = 0.0;
    let mut qual_reads: usize = 0;
    let mut mod_summary = ModSummary::default();
    let mut record = bam::Record::new();

    while let Some(res) = bam.read(&mut record) {
        res.expect("Failure parsing Bam file");
        if (record.flags() & BAM_SECONDARY_OR_SUPPLEMENTARY) != 0 {
            continue;
        }
        let len = record.seq_len() as u32;
        lengths.push(len);
        yield_bases += u64::from(len);
        let aux = record_aux_slice(&record);
        let (aux_err, mask, unknown) = scan_aux_slice(aux);
        if let Some(err) = aux_err.or_else(|| mean_qual_error(record.qual())) {
            qual_sum += err;
            qual_reads += 1;
        }
        mod_summary.merge(mask, unknown);
    }

    lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    Data {
        lengths,
        yield_bases,
        qual_sum,
        qual_reads,
        modifications: mod_summary.format(),
    }
}

#[cfg(unix)]
fn extract_parallel_bam(path: &str, threads: usize, n_targets: i32, voffset: i64) -> Data {
    let file = File::open(Path::new(path)).expect("Error opening BAM file");
    let file_size = file.metadata().expect("Failed reading BAM metadata").len();
    let first_block_offset = (voffset as u64) >> 16;
    let first_in_block_offset = (voffset & 0xFFFF) as usize;

    if first_block_offset >= file_size {
        return Data {
            lengths: Vec::new(),
            yield_bases: 0,
            qual_sum: 0.0,
            qual_reads: 0,
            modifications: "None".to_string(),
        };
    }

    let num_chunks = (threads * 4).max(1);
    let data_len = file_size - first_block_offset;
    let mut chunk_starts = Vec::with_capacity(num_chunks);
    chunk_starts.push(first_block_offset);

    for i in 1..num_chunks {
        let target = first_block_offset + (data_len * i as u64) / num_chunks as u64;
        if let Some(snapped) = find_next_bgzf_block(&file, target, file_size) {
            if snapped < file_size {
                chunk_starts.push(snapped);
            }
        }
    }
    chunk_starts.sort_unstable();
    chunk_starts.dedup();

    let chunk_outputs: Vec<ChunkOutput> = (0..chunk_starts.len())
        .into_par_iter()
        .map(|idx| {
            let chunk_start = chunk_starts[idx];
            let chunk_end = chunk_starts.get(idx + 1).copied().unwrap_or(file_size);
            let initial_offset = if idx == 0 { first_in_block_offset } else { 0 };
            process_bam_chunk(
                &file,
                file_size,
                chunk_start,
                chunk_end,
                initial_offset,
                idx > 0,
                n_targets,
            )
        })
        .collect();

    let total_reads: usize = chunk_outputs.iter().map(|c| c.lengths.len()).sum();
    let mut lengths = Vec::with_capacity(total_reads);
    let mut yield_bases: u64 = 0;
    let mut qual_sum: f64 = 0.0;
    let mut qual_reads: usize = 0;
    let mut mod_summary = ModSummary::default();

    for mut chunk in chunk_outputs {
        lengths.append(&mut chunk.lengths);
        yield_bases += chunk.yield_bases;
        qual_sum += chunk.qual_sum;
        qual_reads += chunk.qual_reads;
        mod_summary.merge_summary(chunk.mods);
    }

    lengths.par_sort_unstable_by(|a, b| b.cmp(a));
    Data {
        lengths,
        yield_bases,
        qual_sum,
        qual_reads,
        modifications: mod_summary.format(),
    }
}

#[cfg(unix)]
fn process_bam_chunk(
    file: &File,
    file_size: u64,
    chunk_start: u64,
    chunk_end: u64,
    initial_offset: usize,
    needs_sync: bool,
    n_targets: i32,
) -> ChunkOutput {
    let mut out = ChunkOutput::default();
    let mut reader = if needs_sync {
        let mut sync = SyncStream::new(file, file_size, chunk_start);
        let synced_pos = loop {
            let p = sync.stream.cursor;
            if !sync.ensure_uncomp_bytes(p + 36) {
                return out;
            }
            match sync.stream.block_file_offset_at(p) {
                Some(off) if off < chunk_end => {}
                _ => return out,
            }
            if validate_record_at(&mut sync, p, n_targets, 0) {
                break p;
            }
            sync.stream.cursor = p + 1;
            if sync.stream.cursor >= STREAM_COMPACT_THRESHOLD {
                sync.stream.compact();
            }
        };
        sync.into_block_reader(synced_pos)
    } else {
        let mut r = BlockReader::new(file, file_size, chunk_start);
        if !r.skip_uncomp_bytes(initial_offset) {
            return out;
        }
        r
    };

    let mut hdr_buf = [0u8; 36];
    let mut mm_scratch = Vec::with_capacity(4096);
    let mut has_rq_or_qs = false;

    loop {
        if !reader.ensure_nonempty_block() {
            break;
        }
        if reader.curr_block_file_off >= chunk_end {
            break;
        }

        let hdr: &[u8] = if reader.avail() >= 36 {
            let s = &reader.block_buf[reader.block_pos..reader.block_pos + 36];
            reader.block_pos += 36;
            s
        } else {
            if !reader.read_exact(&mut hdr_buf) {
                panic!("Truncated BAM record header in chunk starting at {chunk_start}");
            }
            &hdr_buf
        };

        let block_len = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as usize;
        let flag = u16::from_le_bytes([hdr[18], hdr[19]]);
        if (flag & BAM_SECONDARY_OR_SUPPLEMENTARY) != 0 {
            if block_len > 32 && !reader.skip_uncomp_bytes(block_len - 32) {
                panic!("Truncated BAM record body in chunk starting at {chunk_start}");
            }
            continue;
        }

        let l_qname = hdr[12] as usize;
        let n_cigar = u16::from_le_bytes([hdr[16], hdr[17]]) as usize;
        let l_qseq = i32::from_le_bytes([hdr[20], hdr[21], hdr[22], hdr[23]]) as usize;
        let skip_to_qual = l_qname + n_cigar * 4 + ((l_qseq + 1) >> 1);
        let aux_len = (block_len - 32) - skip_to_qual - l_qseq;

        out.lengths.push(l_qseq as u32);
        out.yield_bases += l_qseq as u64;

        if skip_to_qual > 0 && !reader.skip_uncomp_bytes(skip_to_qual) {
            panic!("Truncated BAM record body in chunk starting at {chunk_start}");
        }

        if l_qseq + aux_len <= reader.avail() {
            let qual_start = reader.block_pos;
            let aux_start = qual_start + l_qseq;
            let (aux_err, mask, unknown) = if aux_len > 0 {
                scan_aux_slice(&reader.block_buf[aux_start..aux_start + aux_len])
            } else {
                (None, 0, Vec::new())
            };
            if aux_err.is_some() {
                has_rq_or_qs = true;
            }
            let err_opt = aux_err.or_else(|| mean_qual_error(&reader.block_buf[qual_start..aux_start]));
            if let Some(err) = err_opt {
                out.qual_sum += err;
                out.qual_reads += 1;
            }
            out.mods.merge(mask, unknown);
            reader.block_pos = aux_start + aux_len;
        } else if has_rq_or_qs && aux_len > 0 {
            if l_qseq > 0 && !reader.skip_uncomp_bytes(l_qseq) {
                panic!("Truncated BAM quality array in chunk starting at {chunk_start}");
            }
            let (aux_err, mask, unknown) = scan_aux_streaming(&mut reader, aux_len, &mut mm_scratch);
            if let Some(err) = aux_err {
                out.qual_sum += err;
                out.qual_reads += 1;
            }
            out.mods.merge(mask, unknown);
        } else {
            let qual_err = stream_qual_error(&mut reader, l_qseq);
            let (aux_err, mask, unknown) = if aux_len > 0 {
                scan_aux_streaming(&mut reader, aux_len, &mut mm_scratch)
            } else {
                (None, 0, Vec::new())
            };
            if aux_err.is_some() {
                has_rq_or_qs = true;
            }
            if let Some(err) = aux_err.or(qual_err) {
                out.qual_sum += err;
                out.qual_reads += 1;
            }
            out.mods.merge(mask, unknown);
        }
    }

    out
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

#[cfg(unix)]
fn stream_qual_error(reader: &mut BlockReader<'_>, l_qseq: usize) -> Option<f64> {
    if l_qseq == 0 || !reader.ensure_nonempty_block() {
        return None;
    }
    if reader.block_buf[reader.block_pos] == 0xFF {
        reader.skip_uncomp_bytes(l_qseq);
        return None;
    }
    let table = &*PHRED_TO_ERROR;
    let mut rem = l_qseq;
    let mut err_sum = 0.0_f64;
    while rem > 0 {
        if !reader.ensure_nonempty_block() {
            return None;
        }
        let take = rem.min(reader.avail());
        let slice = &reader.block_buf[reader.block_pos..reader.block_pos + take];
        let mut chunk_sum = 0.0_f64;
        for &q in slice {
            chunk_sum += table[q as usize];
        }
        err_sum += chunk_sum;
        reader.block_pos += take;
        rem -= take;
    }
    Some(err_sum / l_qseq as f64)
}

#[inline]
fn aux_fixed_size(val_type: u8) -> Option<usize> {
    match val_type {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        b'd' => Some(8),
        _ => None,
    }
}

#[inline]
fn parse_aux_number(val_type: u8, payload: &[u8]) -> Option<f64> {
    match val_type {
        b'c' => Some(f64::from(payload[0] as i8)),
        b'C' => Some(f64::from(payload[0])),
        b's' => Some(f64::from(i16::from_le_bytes([payload[0], payload[1]]))),
        b'S' => Some(f64::from(u16::from_le_bytes([payload[0], payload[1]]))),
        b'i' => Some(f64::from(i32::from_le_bytes([
            payload[0], payload[1], payload[2], payload[3],
        ]))),
        b'I' => Some(f64::from(u32::from_le_bytes([
            payload[0], payload[1], payload[2], payload[3],
        ]))),
        b'f' => Some(f64::from(f32::from_le_bytes([
            payload[0], payload[1], payload[2], payload[3],
        ]))),
        b'd' => Some(f64::from_le_bytes([
            payload[0], payload[1], payload[2], payload[3], payload[4], payload[5], payload[6],
            payload[7],
        ])),
        _ => None,
    }
}

#[inline]
fn update_rq_qs_err(
    tag0: u8,
    tag1: u8,
    val_type: u8,
    payload: &[u8],
    rq_err: &mut Option<f64>,
    qs_err: &mut Option<f64>,
) {
    if tag0 == b'r' && tag1 == b'q' && matches!(val_type, b'f' | b'd') {
        if let Some(v) = parse_aux_number(val_type, payload) {
            if v.is_finite() && (0.0..=1.0).contains(&v) {
                *rq_err = Some(1.0 - v);
            }
        }
    } else if tag0 == b'q' && tag1 == b's' && val_type != b'A' {
        if let Some(v) = parse_aux_number(val_type, payload) {
            if v.is_finite() && v >= 0.0 {
                *qs_err = Some(10.0_f64.powf(-v / 10.0));
            }
        }
    }
}

#[cfg(unix)]
fn scan_aux_streaming(
    reader: &mut BlockReader<'_>,
    aux_len: usize,
    mm_scratch: &mut Vec<u8>,
) -> (Option<f64>, u32, Vec<String>) {
    let mut rq_err: Option<f64> = None;
    let mut qs_err: Option<f64> = None;
    let mut found_upper_mm = false;
    let mut mod_result = (0u32, Vec::new());

    let mut rem = aux_len;
    let mut tag_hdr = [0u8; 3];
    let mut arr_hdr = [0u8; 5];
    let mut num_buf = [0u8; 8];

    while rem >= 3 {
        if !reader.read_exact(&mut tag_hdr) {
            break;
        }
        rem -= 3;
        let tag0 = tag_hdr[0];
        let tag1 = tag_hdr[1];
        let val_type = tag_hdr[2];

        if let Some(size) = aux_fixed_size(val_type) {
            if rem < size {
                break;
            }
            let is_rq_or_qs = (tag0 == b'r' && tag1 == b'q') || (tag0 == b'q' && tag1 == b's');
            if is_rq_or_qs && val_type != b'A' {
                if !reader.read_exact(&mut num_buf[..size]) {
                    break;
                }
                rem -= size;
                update_rq_qs_err(
                    tag0,
                    tag1,
                    val_type,
                    &num_buf[..size],
                    &mut rq_err,
                    &mut qs_err,
                );
                if rq_err.is_some() && found_upper_mm {
                    if rem > 0 {
                        reader.skip_uncomp_bytes(rem);
                    }
                    return (rq_err, mod_result.0, mod_result.1);
                }
            } else {
                if !reader.skip_uncomp_bytes(size) {
                    break;
                }
                rem -= size;
            }
        } else {
            match val_type {
                b'Z' | b'H' => {
                    if val_type == b'Z'
                        && ((tag0 == b'M' && tag1 == b'M')
                            || (tag0 == b'M' && tag1 == b'm' && !found_upper_mm))
                    {
                        let consumed = reader.read_z_into(mm_scratch, rem);
                        rem -= consumed;
                        mod_result = parse_mm_bytes(mm_scratch);
                        if tag1 == b'M' {
                            found_upper_mm = true;
                            if rq_err.is_some() {
                                if rem > 0 {
                                    reader.skip_uncomp_bytes(rem);
                                }
                                return (rq_err, mod_result.0, mod_result.1);
                            }
                        }
                    } else {
                        let consumed = reader.skip_until_nul(rem);
                        rem -= consumed;
                    }
                }
                b'B' => {
                    if rem < 5 || !reader.read_exact(&mut arr_hdr) {
                        break;
                    }
                    rem -= 5;
                    let Some(elem_size) = (match arr_hdr[0] {
                        b'c' | b'C' => Some(1usize),
                        b's' | b'S' => Some(2usize),
                        b'i' | b'I' | b'f' => Some(4usize),
                        _ => None,
                    }) else {
                        break;
                    };
                    let count = u32::from_le_bytes([arr_hdr[1], arr_hdr[2], arr_hdr[3], arr_hdr[4]])
                        as usize;
                    let Some(byte_len) = count.checked_mul(elem_size) else {
                        break;
                    };
                    if byte_len > rem || !reader.skip_uncomp_bytes(byte_len) {
                        break;
                    }
                    rem -= byte_len;
                }
                _ => break,
            }
        }
    }

    if rem > 0 {
        reader.skip_uncomp_bytes(rem);
    }

    (rq_err.or(qs_err), mod_result.0, mod_result.1)
}

#[cfg(unix)]
fn read_exact_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> std::io::Result<()> {
    while !buf.is_empty() {
        let n = file.read_at(buf, offset)?;
        if n == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Unexpected EOF",
            ));
        }
        buf = &mut buf[n..];
        offset += n as u64;
    }
    Ok(())
}

/// Parses a BGZF header returning `(bsize, header_len)` where `header_len = 12 + xlen`.
#[inline]
fn parse_bgzf_header(h: &[u8]) -> Option<(usize, usize)> {
    if h.len() < 18 {
        return None;
    }
    if h[0] != 0x1f || h[1] != 0x8b || h[2] != 0x08 || (h[3] & 0x04) == 0 {
        return None;
    }
    let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
    if xlen != 6 {
        return None;
    }
    if h[12] != b'B' || h[13] != b'C' || u16::from_le_bytes([h[14], h[15]]) != 2 {
        return None;
    }
    let bsize = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
    let header_len = 12 + xlen;
    if bsize < header_len + 8 {
        return None;
    }
    Some((bsize, header_len))
}

#[cfg(unix)]
fn find_next_bgzf_block(file: &File, start_off: u64, file_size: u64) -> Option<u64> {
    if start_off + 26 > file_size {
        return None;
    }
    const WIN_SIZE: usize = 128 * 1024;
    let mut win = vec![0u8; WIN_SIZE];
    let mut off = start_off;

    while off + 18 <= file_size {
        let to_read = ((file_size - off).min(WIN_SIZE as u64)) as usize;
        if read_exact_at(file, &mut win[..to_read], off).is_err() {
            return None;
        }
        let search_end = to_read.saturating_sub(17);
        let mut i = 0;
        while i < search_end {
            let Some(rel) = win[i..search_end].iter().position(|&b| b == 0x1f) else {
                break;
            };
            i += rel;
            if let Some((bsize, _)) = parse_bgzf_header(&win[i..to_read]) {
                let cand = off + i as u64;
                let next_off = cand + bsize as u64;
                if next_off == file_size {
                    return Some(cand);
                }
                if next_off + 18 <= file_size {
                    let next_in_win = i + bsize;
                    let valid_next = if next_in_win + 18 <= to_read {
                        parse_bgzf_header(&win[next_in_win..next_in_win + 18])
                            .is_some_and(|(nb, _)| next_off + nb as u64 <= file_size)
                    } else {
                        let mut nh = [0u8; 18];
                        read_exact_at(file, &mut nh, next_off).is_ok()
                            && parse_bgzf_header(&nh)
                                .is_some_and(|(nb, _)| next_off + nb as u64 <= file_size)
                    };
                    if valid_next {
                        return Some(cand);
                    }
                }
            }
            i += 1;
        }
        if off + to_read as u64 >= file_size {
            break;
        }
        off += search_end as u64;
    }

    None
}

#[cfg(unix)]
struct RecordHeader {
    block_len: usize,
    l_qname: usize,
    n_cigar: usize,
    aux_rel: usize,
}

#[cfg(unix)]
#[inline]
fn parse_record_header(hdr: &[u8], n_targets: i32) -> Option<RecordHeader> {
    if hdr.len() < 36 {
        return None;
    }
    let block_len = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as usize;
    if !(32..=MAX_BAM_RECORD_LEN).contains(&block_len) {
        return None;
    }
    let tid = i32::from_le_bytes([hdr[4], hdr[5], hdr[6], hdr[7]]);
    if tid < -1 || tid >= n_targets {
        return None;
    }
    let pos = i32::from_le_bytes([hdr[8], hdr[9], hdr[10], hdr[11]]);
    if pos < -1 {
        return None;
    }
    let l_qname = hdr[12] as usize;
    if l_qname < 2 {
        return None;
    }
    let n_cigar = u16::from_le_bytes([hdr[16], hdr[17]]) as usize;
    let flag = u16::from_le_bytes([hdr[18], hdr[19]]);
    if flag > 0x0FFF {
        return None;
    }
    let l_qseq_i32 = i32::from_le_bytes([hdr[20], hdr[21], hdr[22], hdr[23]]);
    if !(0..=MAX_BAM_SEQ_LEN).contains(&l_qseq_i32) {
        return None;
    }
    let l_qseq = l_qseq_i32 as usize;
    let mtid = i32::from_le_bytes([hdr[24], hdr[25], hdr[26], hdr[27]]);
    if mtid < -1 || mtid >= n_targets {
        return None;
    }
    let mpos = i32::from_le_bytes([hdr[28], hdr[29], hdr[30], hdr[31]]);
    if mpos < -1 {
        return None;
    }
    let seq_bytes = (l_qseq + 1) >> 1;
    let fixed_payload = l_qname
        .checked_add(n_cigar * 4)?
        .checked_add(seq_bytes)?
        .checked_add(l_qseq)?;
    if 32 + fixed_payload > block_len {
        return None;
    }
    Some(RecordHeader {
        block_len,
        l_qname,
        n_cigar,
        aux_rel: 36 + fixed_payload,
    })
}

#[cfg(unix)]
fn validate_aux_framing(aux: &[u8]) -> bool {
    let mut i = 0;
    while i < aux.len() {
        if i + 3 > aux.len() {
            return false;
        }
        if !aux[i].is_ascii_alphabetic() || !aux[i + 1].is_ascii_alphanumeric() {
            return false;
        }
        let payload_start = i + 3;
        let next = if let Some(size) = aux_fixed_size(aux[i + 2]) {
            payload_start + size
        } else {
            match aux[i + 2] {
                b'Z' | b'H' => {
                    let Some(nul_rel) = aux[payload_start..].iter().position(|&b| b == 0) else {
                        return false;
                    };
                    payload_start + nul_rel + 1
                }
                b'B' => {
                    if payload_start + 5 > aux.len() {
                        return false;
                    }
                    let elem_size = match aux[payload_start] {
                        b'c' | b'C' => 1usize,
                        b's' | b'S' => 2usize,
                        b'i' | b'I' | b'f' => 4usize,
                        _ => return false,
                    };
                    let count = u32::from_le_bytes([
                        aux[payload_start + 1],
                        aux[payload_start + 2],
                        aux[payload_start + 3],
                        aux[payload_start + 4],
                    ]) as usize;
                    let Some(byte_len) = count.checked_mul(elem_size) else {
                        return false;
                    };
                    let Some(end) = (payload_start + 5).checked_add(byte_len) else {
                        return false;
                    };
                    end
                }
                _ => return false,
            }
        };
        if next > aux.len() {
            return false;
        }
        i = next;
    }
    i == aux.len()
}

#[cfg(unix)]
fn validate_record_at(sync: &mut SyncStream<'_>, p: usize, n_targets: i32, depth: u8) -> bool {
    if !sync.ensure_uncomp_bytes(p + 36) {
        return false;
    }
    let Some(hdr) = parse_record_header(&sync.stream.data[p..p + 36], n_targets) else {
        return false;
    };

    // Validate QNAME before pulling large records into the stream buffer.
    let qname_end = p + 36 + hdr.l_qname;
    if !sync.ensure_uncomp_bytes(qname_end) {
        return false;
    }
    let qname = &sync.stream.data[p + 36..qname_end];
    if qname[hdr.l_qname - 1] != 0
        || !qname[..hdr.l_qname - 1]
            .iter()
            .all(|&b| (0x21..=0x7e).contains(&b))
    {
        return false;
    }

    let rec_end = p + 4 + hdr.block_len;
    if !sync.ensure_uncomp_bytes(rec_end) {
        return false;
    }

    let cigar_start = qname_end;
    let cigar_end = cigar_start + hdr.n_cigar * 4;
    for op_bytes in sync.stream.data[cigar_start..cigar_end].chunks_exact(4) {
        let raw = u32::from_le_bytes([op_bytes[0], op_bytes[1], op_bytes[2], op_bytes[3]]);
        if (raw & 0x0f) > 8 || (raw >> 4) == 0 {
            return false;
        }
    }

    if !validate_aux_framing(&sync.stream.data[p + hdr.aux_rel..rec_end]) {
        return false;
    }

    if depth == 0 {
        if !sync.ensure_uncomp_bytes(rec_end + 1) {
            return true;
        }
        return validate_record_at(sync, rec_end, n_targets, 1);
    }

    true
}

#[inline]
fn record_aux_slice(record: &bam::Record) -> &[u8] {
    let inner = record.inner();
    let core = &inner.core;
    let aux_offset = ((core.n_cigar as usize) << 2)
        + (core.l_qname as usize)
        + (((core.l_qseq as usize) + 1) >> 1)
        + (core.l_qseq as usize);
    let l_data = inner.l_data as usize;
    if inner.data.is_null() || aux_offset >= l_data {
        return &[];
    }
    unsafe { std::slice::from_raw_parts(inner.data.add(aux_offset), l_data - aux_offset) }
}

/// Scans BAM raw auxiliary bytes in a single pass for `rq`, `qs`, and `MM`/`Mm` tags.
fn scan_aux_slice(aux: &[u8]) -> (Option<f64>, u32, Vec<String>) {
    let mut rq_err: Option<f64> = None;
    let mut qs_err: Option<f64> = None;
    let mut found_upper_mm = false;
    let mut mod_result = (0u32, Vec::new());

    let mut i = 0;
    while i + 3 <= aux.len() {
        let tag0 = aux[i];
        let tag1 = aux[i + 1];
        let val_type = aux[i + 2];
        let payload_start = i + 3;

        let next = if let Some(size) = aux_fixed_size(val_type) {
            let end = payload_start + size;
            if end > aux.len() {
                break;
            }
            if val_type != b'A' {
                update_rq_qs_err(
                    tag0,
                    tag1,
                    val_type,
                    &aux[payload_start..end],
                    &mut rq_err,
                    &mut qs_err,
                );
                if rq_err.is_some() && found_upper_mm {
                    break;
                }
            }
            end
        } else {
            match val_type {
                b'Z' | b'H' => {
                    let Some(nul_rel) = aux[payload_start..].iter().position(|&b| b == 0) else {
                        break;
                    };
                    let str_end = payload_start + nul_rel;
                    if val_type == b'Z' {
                        if tag0 == b'M' && tag1 == b'M' {
                            mod_result = parse_mm_bytes(&aux[payload_start..str_end]);
                            found_upper_mm = true;
                            if rq_err.is_some() {
                                break;
                            }
                        } else if tag0 == b'M' && tag1 == b'm' && !found_upper_mm {
                            mod_result = parse_mm_bytes(&aux[payload_start..str_end]);
                        }
                    }
                    str_end + 1
                }
                b'B' => {
                    if payload_start + 5 > aux.len() {
                        break;
                    }
                    let elem_size = match aux[payload_start] {
                        b'c' | b'C' => 1usize,
                        b's' | b'S' => 2usize,
                        b'i' | b'I' | b'f' => 4usize,
                        _ => break,
                    };
                    let count = u32::from_le_bytes([
                        aux[payload_start + 1],
                        aux[payload_start + 2],
                        aux[payload_start + 3],
                        aux[payload_start + 4],
                    ]) as usize;
                    let Some(byte_len) = count.checked_mul(elem_size) else {
                        break;
                    };
                    let Some(end) = (payload_start + 5).checked_add(byte_len) else {
                        break;
                    };
                    end
                }
                _ => break,
            }
        };

        if next > aux.len() {
            break;
        }
        i = next;
    }

    (rq_err.or(qs_err), mod_result.0, mod_result.1)
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
    use rust_htslib::bam::record::Aux;
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
    fn scan_aux_handles_rq_qs_priority_and_qual_fallback() {
        let mut rec_mm = make_record(b"ACGTACGT", Some(0.99), Some((b"MM", "C+mh?,0;")));
        rec_mm
            .push_aux(b"ML", Aux::ArrayU8((&[250u8, 12u8][..]).into()))
            .unwrap();
        let (err_opt, mask, unknown) = scan_aux_slice(record_aux_slice(&rec_mm));
        assert!((err_opt.unwrap() - 0.01).abs() < 1e-6);
        assert_eq!((mask, unknown), (MOD_5MC | MOD_5HMC, Vec::new()));

        let mut rec_qs_float = make_record(b"ACGTACGT", None, Some((b"Mm", "A+a.,0;")));
        rec_qs_float.push_aux(b"qs", Aux::Float(20.0)).unwrap();
        let (err_opt, mask, unknown) = scan_aux_slice(record_aux_slice(&rec_qs_float));
        assert!((err_opt.unwrap() - 0.01).abs() < 1e-6);
        assert_eq!((mask, unknown), (MOD_6MA, Vec::new()));

        let mut rec_qs_then_rq = make_record(b"ACGTACGT", None, None);
        rec_qs_then_rq.push_aux(b"qs", Aux::I32(10)).unwrap();
        rec_qs_then_rq.push_aux(b"rq", Aux::Float(0.99)).unwrap();
        let (err_opt, _, _) = scan_aux_slice(record_aux_slice(&rec_qs_then_rq));
        assert!((err_opt.unwrap() - 0.01).abs() < 1e-6);

        let rec_none = make_record(b"ACGTACGT", None, None);
        assert_eq!(
            scan_aux_slice(record_aux_slice(&rec_none)),
            (None, 0, Vec::new())
        );
        assert!((mean_qual_error(rec_none.qual()).unwrap() - 0.001).abs() < 1e-9);
        assert_eq!(mean_qual_error(&[0xFF, 0xFF, 0xFF]), None);
    }

    #[test]
    fn extract_parallel_matches_streaming_and_filters_secondary_supplementary() {
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

            // Illumina-like read without rq or qs, falling back to base QUAL (Q20 -> 0.01).
            let mut r3 = bam::Record::new();
            r3.set(b"illumina_read", None, b"ACGT", &[20, 20, 20, 20]);

            // Read with negative rq (-1.0) and missing QUAL (*): ignored in qual_reads.
            let mut r4 = bam::Record::new();
            r4.set(b"unscored", None, b"ACGTAC", &[0xFF; 6]);
            r4.push_aux(b"rq", Aux::Float(-1.0)).unwrap();

            // Secondary (0x100) and supplementary (0x800) alignments should be filtered out.
            let mut r_sec = make_record(b"TTTTTTTTTTTTTTTTTTTT", Some(0.5), Some((b"MM", "G+o?,0;")));
            r_sec.set_flags(0x100);
            let mut r_sup = make_record(b"CCCCCCCCCCCCCCCCCCCC", Some(0.5), Some((b"MM", "T+g?,0;")));
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
        let stream_reader = bam::Reader::from_path(&cli.input).unwrap();
        let stream_data = extract_streaming(&cli, stream_reader);
        std::fs::remove_file(&bam_path).ok();

        assert_eq!(data.lengths, stream_data.lengths);
        assert_eq!(data.yield_bases, stream_data.yield_bases);
        assert_eq!(data.qual_reads, stream_data.qual_reads);
        assert!((data.qual_sum - stream_data.qual_sum).abs() < 1e-9);
        assert_eq!(data.modifications, stream_data.modifications);
        assert_eq!(data.lengths, vec![12, 10, 6, 4]);
        assert_eq!(data.yield_bases, 32);
        assert_eq!(data.qual_reads, 3);
        assert!((data.qual_sum - (0.01 + 0.001 + 0.01)).abs() < 1e-6);
        assert_eq!(data.modifications, "5mC;5hmC;4mC;6mA");
        assert_eq!(parsed_header.to_hashmap()["RG"][0]["PU"], "flowcell_pu");

        for repo_bam in ["test-data/small-test-ubam.bam", "test-data/small-test-phased.bam"] {
            let repo_cli = crate::Cli {
                input: repo_bam.to_string(),
                threads: 4,
                reference: None,
            };
            let (par_data, _) = extract(&repo_cli);
            let stream_reader = bam::Reader::from_path(repo_bam).unwrap();
            let seq_data = extract_streaming(&repo_cli, stream_reader);
            assert_eq!(par_data.lengths, seq_data.lengths);
            assert_eq!(par_data.yield_bases, seq_data.yield_bases);
            assert_eq!(par_data.qual_reads, seq_data.qual_reads);
            assert!((par_data.qual_sum - seq_data.qual_sum).abs() < 1e-9);
            assert_eq!(par_data.modifications, seq_data.modifications);
        }
    }

    #[test]
    fn extract_detects_modifications_across_many_bgzf_blocks_and_chunks() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let bam_path = std::env::temp_dir().join(format!("cramino_multiblock_{nanos}.bam"));

        let header = Header::new();
        let total_reads = 6000;
        {
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut plain = bam::Record::new();
            let mut pseudo_seq = vec![0u8; 512];
            for (idx, base) in pseudo_seq.iter_mut().enumerate() {
                *base = b"ACGT"[(idx * 7 + 3) % 4];
            }
            let qual = vec![30u8; 512];
            plain.set(b"plain_long", None, &pseudo_seq, &qual);
            plain.push_aux(b"rq", Aux::Float(0.99)).unwrap();

            let with_5mc = make_record(&pseudo_seq[..250], Some(0.99), Some((b"MM", "C+m?,0;")));
            let with_6ma_tail = make_record(&pseudo_seq[..300], Some(0.99), Some((b"MM", "A+a.,0;")));

            writer.write(&with_5mc).unwrap();
            for i in 1..(total_reads - 1) {
                let name = format!("read_{i:05}");
                plain.set_qname(name.as_bytes());
                writer.write(&plain).unwrap();
            }
            writer.write(&with_6ma_tail).unwrap();
        }

        let cli = crate::Cli {
            input: bam_path.to_str().unwrap().to_string(),
            threads: 8,
            reference: None,
        };
        let (par_data, _) = extract(&cli);
        let stream_reader = bam::Reader::from_path(&cli.input).unwrap();
        let seq_data = extract_streaming(&cli, stream_reader);
        std::fs::remove_file(&bam_path).ok();

        assert_eq!(par_data.lengths, seq_data.lengths);
        assert_eq!(par_data.yield_bases, seq_data.yield_bases);
        assert_eq!(par_data.qual_reads, seq_data.qual_reads);
        assert!((par_data.qual_sum - seq_data.qual_sum).abs() < 1e-6);
        assert_eq!(par_data.lengths.len(), total_reads);
        assert_eq!(par_data.modifications, "5mC;6mA");
        assert_eq!(par_data.modifications, seq_data.modifications);
    }

    #[test]
    fn extract_handles_records_and_large_ml_arrays_spanning_blocks() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let bam_path = std::env::temp_dir().join(format!("cramino_longread_{nanos}.bam"));

        let header = Header::new();
        let mut long_seq = vec![0u8; 150_000];
        let mut state: u32 = 0x1234_5678;
        for base in &mut long_seq {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            *base = b"ACGT"[((state >> 24) & 3) as usize];
        }
        let long_qual = vec![20u8; long_seq.len()];
        let large_ml: Vec<u8> = (0..120_000).map(|i| (i % 256) as u8).collect();

        {
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let mut r_short = make_record(b"ACGTACGT", Some(0.99), None);
            let mut r_long = bam::Record::new();
            r_long.set(b"ultra_long_read", None, &long_seq, &long_qual);
            r_long
                .push_aux(b"ML", Aux::ArrayU8((&large_ml[..]).into()))
                .unwrap();
            r_long.push_aux(b"MM", Aux::String("C+h?,0;")).unwrap();
            r_long.push_aux(b"rq", Aux::Float(0.995)).unwrap();

            // Also include a cross-block read with base QUAL fallback only (no rq/qs).
            let mut r_qual_only = bam::Record::new();
            r_qual_only.set(b"qual_only_long", None, &long_seq, &long_qual);

            for i in 0..10 {
                r_short.set_qname(format!("pre_{i}").as_bytes());
                writer.write(&r_short).unwrap();
            }
            writer.write(&r_long).unwrap();
            writer.write(&r_qual_only).unwrap();
            for i in 0..10 {
                r_short.set_qname(format!("post_{i}").as_bytes());
                writer.write(&r_short).unwrap();
            }
        }

        let cli = crate::Cli {
            input: bam_path.to_str().unwrap().to_string(),
            threads: 8,
            reference: None,
        };
        let (par_data, _) = extract(&cli);
        let stream_reader = bam::Reader::from_path(&cli.input).unwrap();
        let seq_data = extract_streaming(&cli, stream_reader);
        std::fs::remove_file(&bam_path).ok();

        assert_eq!(par_data.lengths, seq_data.lengths);
        assert_eq!(par_data.yield_bases, seq_data.yield_bases);
        assert_eq!(par_data.qual_reads, seq_data.qual_reads);
        assert!((par_data.qual_sum - seq_data.qual_sum).abs() < 1e-6);
        assert_eq!(par_data.modifications, "5hmC");
        assert_eq!(par_data.lengths[0], 150_000);
        assert_eq!(par_data.lengths.len(), 22);
    }
}
