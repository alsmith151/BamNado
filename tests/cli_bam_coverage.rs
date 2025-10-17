use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn test_data_dir() -> PathBuf {
    workspace_root().join("test/data")
}

#[test]
fn bam_coverage_writes_bedgraph_by_default() {
    let bam = test_data_dir().join("test.bam");
    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");

    // default output is input.bedgraph; use a temp dir to avoid polluting repo
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("out.bedgraph");

    cmd.arg("bam-coverage")
        .arg("--bam").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--bin-size").arg("1000")
        .arg("--norm-method").arg("raw")
        .arg("--use-fragment");

    cmd.assert()
        .success()
        .stderr(predicate::str::contains("Successfully wrote output"));

    // verify file exists and has some content
    assert!(out.exists(), "bedgraph output missing");
    let meta = std::fs::metadata(&out).expect("metadata");
    assert!(meta.len() > 0, "bedgraph is empty");
}

#[test]
fn bam_coverage_writes_bigwig() {
    let bam = test_data_dir().join("test.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("signal.bw");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("bam-coverage")
        .arg("--bam").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--bin-size").arg("1000")
        .arg("--norm-method").arg("raw");

    cmd.assert().success();
    assert!(out.exists(), "bigwig output missing");
}
