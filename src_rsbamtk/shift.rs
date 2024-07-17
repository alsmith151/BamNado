//use rust_htslib::bam::record::Cigar;
use log::{info, warn};
use rust_htslib::bam::{Format, Header, Read};
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;

// Copying from this:
// def shiftRead(b, chromDict, args):
//     if not b.is_proper_pair:
//         return None
//     tLen = getTLen(b, notAbs=True)
//     start = b.pos
//     end = start + b.query_alignment_end
//     if b.is_reverse and not b.is_read2:
//         end -= args.shift[2]
//         deltaTLen = args.shift[3] - args.shift[2]
//     elif b.is_reverse and b.is_read2:
//         end += args.shift[1]
//         deltaTLen = args.shift[1] - args.shift[0]
//     elif not b.is_reverse and not b.is_read2:
//         start += args.shift[0]
//         deltaTLen = args.shift[1] - args.shift[0]
//     else:
//         start -= args.shift[3]
//         deltaTLen = args.shift[3] - args.shift[2]

//     # Sanity check
//     if end - start < 1:
//         if b.is_reverse:
//             start = end - 1
//         else:
//             end = start + 1
//     if start < 0:
//         start = 0
//     if end > chromDict[b.reference_name]:
//         end = chromDict[b.reference_name]
//     if end - start < 1:
//         return None

//     # create a new read
//     b2 = pysam.AlignedSegment()
//     b2.query_name = b.query_name
//     b2.flag = b.flag
//     b2.reference_id = b.reference_id
//     b2.reference_start = start
//     b2.mapping_quality = b.mapping_quality
//     b2.cigar = ((0, end - start),)  # Returned cigar is only matches
//     if tLen < 0:
//         b2.template_length = tLen - deltaTLen
//     else:
//         b2.template_length = tLen + deltaTLen
//     b2.next_reference_id = b.next_reference_id
//     b2.next_reference_start = b.next_reference_start
//     if b.is_proper_pair:
//         if b2.is_read2 and b2.is_reverse:
//             b2.next_reference_start += args.shift[0]
//         elif not b2.is_read2 and b2.is_reverse:
//             b2.next_reference_start -= args.shift[3]

//     return b2

//# Sanity check
//     if end - start < 1:
//         if b.is_reverse:
//             start = end - 1
//         else:
//             end = start + 1
//     if start < 0:
//         start = 0
//     if end > chromDict[b.reference_name]:
//         end = chromDict[b.reference_name]
//     if end - start < 1:
//         return None

fn sanity_check_coordinates(
    mut start: i64,
    mut end: i64,
    reverse: bool,
    chromsize: i64,
) -> Option<(i64, i64)> {
    if end - start < 1 {
        match reverse {
            true => start = end - 1,
            false => end = start + 1,
        }
    }

    if start < 0 {
        start = 0
    }

    if end > chromsize {
        end = chromsize
    }

    if end - start >= 1 {
        Some((start, end))
    } else {
        None
    }
}

fn set_up_chromsizes(
    header: &rust_htslib::bam::HeaderView,
) -> Result<HashMap<u32, u64>, Box<dyn Error>> {
    let tids: HashMap<u32, u64> = header
        .target_names()
        .iter()
        .map(|n| header.tid(n))
        .filter(|tid| tid.is_some())
        .map(|tid| (tid.unwrap(), header.target_len(tid.unwrap()).unwrap()))
        .collect();
    Ok(tids)
}

pub fn atac_shift_bam<P>(bam_input: P, bam_output: P) -> Result<(), rust_htslib::errors::Error>
where
    P: AsRef<Path>,
{
    let mut reader = rust_htslib::bam::Reader::from_path(bam_input)?;
    let header = Header::from_template(reader.header());
    let mut writer = rust_htslib::bam::Writer::from_path(bam_output, &header, Format::Bam)?;
    let chrom_dict = set_up_chromsizes(reader.header()).expect("Couldn't read chromsizes");
    let shift = vec![4, -5, 5, -4];

    let mut read_counter = 0;
    for result in reader.records() {
        let mut record = result?;

        if record.is_proper_pair() {
            let mut tlen = record.insert_size();
            let mut start = record.pos();
            let mut end = (start as usize + record.seq_len()) as i64;
            let reverse = record.is_reverse();
            let first_in_template = record.is_first_in_template();
            let chromsize = chrom_dict
                .get(&(record.tid() as u32))
                .expect("Missing chromsize");

            let dtlen = match (reverse, first_in_template) {
                (true, true) => {
                    end += shift[1];
                    shift[1] - shift[0]
                }
                (true, false) => {
                    end -= shift[2];
                    shift[3] - shift[2]
                }
                (false, true) => {
                    start -= shift[3];
                    shift[3] - shift[2]
                }
                (false, false) => {
                    start += shift[0];
                    shift[1] - shift[0]
                }
            };

            if let Some((start, _end)) =
                sanity_check_coordinates(start, end, reverse, *chromsize as i64)
            {
                // Edit the record
                record.set_pos(start);

                if tlen > 0 {
                    tlen += dtlen;
                } else {
                    tlen -= dtlen;
                }
                record.set_insert_size(tlen);

                match (reverse, first_in_template) {
                    (true, true) => {
                        let mpos = record.mpos() + shift[0];
                        record.set_mpos(mpos)
                    }
                    (true, false) => {
                        let mpos = record.mpos() - shift[3];
                        record.set_mpos(mpos)
                    }
                    _ => {}
                };

                writer.write(&record)?;
            }
            // Update counter
            read_counter += 1;
            if read_counter % 100000 == 0 {
                println!("Shifted {} reads", read_counter);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use tempdir::TempDir;

    use crate::shift;

    #[test]
    fn shift_bam() {
        let bam = "test/test.bam";

        {
            let tmp = TempDir::new("shift_bam_test").expect("Failed to make tmpdir");

            let out = tmp.path().to_owned().join("test.bam");

            let result = shift::atac_shift_bam(
                bam,
                out.as_path().to_str().expect("Cannot convert"),
            );
            let out_path = out.exists();
            assert_eq!(result.is_ok(), true);
            assert_eq!(out_path, true);
        }
    }
}
