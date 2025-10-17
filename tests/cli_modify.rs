use assert_cmd::prelude::*;
use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")) }
fn test_data_dir() -> PathBuf { workspace_root().join("test/data") }

#[test]
fn modify_runs_and_writes_bam() {
    let bam = test_data_dir().join("test.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let out = temp.path().join("modified.bam");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("modify")
        .arg("--input").arg(&bam)
        .arg("--output").arg(&out)
        .arg("--tn5-shift");

    cmd.assert().success();
    assert!(out.exists(), "modified BAM missing");
}
