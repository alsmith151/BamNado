use assert_cmd::prelude::*;
use serde::de;
use std::path::PathBuf;
use std::process::Command;
use polars::prelude::*;

fn workspace_root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")) }
fn test_data_dir() -> PathBuf { workspace_root().join("test/data") }

#[test]
fn fragments_creates_gz_output() {
    let bam = test_data_dir().join("test_sc.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("frags.tsv.gz");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("fragments")
        .arg("--input").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--min-fragment-length").arg("50")
        .arg("--max-fragment-length").arg("1000")
        .arg("--chunk-size").arg("5000");

    cmd.assert().success();

    assert!(out.exists(), "fragments output missing");
    let meta = std::fs::metadata(&out).expect("metadata");
    assert!(meta.len() > 0, "fragments output empty");
}


#[test]
fn test_collapse_within_barcode() {
    let bam = test_data_dir().join("test_sc.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("collapsed_frags.tsv.gz");
    let example_path = test_data_dir().join("frags_within.tsv");


    let example = CsvReadOptions::default()
        .with_has_header(false)
        .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
        .try_into_reader_with_file_path(Some(example_path))
        .expect("Failed to read example file")
        .finish()
        .expect("Error creating output frame");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("fragments")
        .arg("--input").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--min-fragment-length").arg("50")
        .arg("--max-fragment-length").arg("1000")
        .arg("--collapse-within-barcode")
        .arg("--chunk-size").arg("5000");

    cmd.assert().success();

    assert!(out.exists(), "collapsed fragments output missing");
    let meta = std::fs::metadata(&out).expect("metadata");
    assert!(meta.len() > 0, "collapsed fragments output empty");

    let df = CsvReadOptions::default()
        .with_has_header(false)
        .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
        .try_into_reader_with_file_path(Some(out))
        .expect("Failed to read output file")
        .finish()
        .expect("Error creating output frame");


    // Check that the output matches expected collapsed fragments
    assert_eq!(df.height(), example.height(), "Number of rows should match");

}