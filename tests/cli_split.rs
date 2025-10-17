use assert_cmd::prelude::*;
use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")) }
fn test_data_dir() -> PathBuf { workspace_root().join("test/data") }

#[test]
fn split_writes_filtered_bam() {
    let bam = test_data_dir().join("test.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("filtered.bam");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("split")
        .arg("--input").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--min-mapq").arg("0")
        .arg("--min-length").arg("20")
        .arg("--max-length").arg("1000");

    cmd.assert().success();
    assert!(out.exists(), "filtered BAM missing");
    let meta = std::fs::metadata(&out).expect("metadata");
    assert!(meta.len() > 0, "filtered BAM empty");
}
