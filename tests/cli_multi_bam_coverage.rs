use assert_cmd::prelude::*;
use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")) }
fn test_data_dir() -> PathBuf { workspace_root().join("test/data") }

#[test]
fn multi_bam_coverage_writes_tsv() {
    let bam = test_data_dir().join("test.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("multi.tsv");
    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");

    cmd.arg("multi-bam-coverage")
        .arg("--bams").arg(&bam).arg(&bam) // use same file twice for simplicity
        .arg("--output").arg(&out)
        .arg("--bin-size").arg("1000")
        .arg("--norm-method").arg("raw");

    cmd.assert().success();

    assert!(out.exists(), "tsv output missing");
    let content = std::fs::read_to_string(&out).expect("read tsv");
    assert!(content.lines().next().unwrap().contains("chrom"), "header present");
}
