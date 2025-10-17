use assert_cmd::prelude::*;
use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")) }
fn test_data_dir() -> PathBuf { workspace_root().join("test/data") }

#[test]
fn split_exogenous_produces_outputs() {
    let bam = test_data_dir().join("test.bam");
    let temp = assert_fs::TempDir::new().unwrap();
    let prefix = temp.path().join("split");

    let mut cmd = Command::cargo_bin("bamnado").expect("binary exists");
    cmd.arg("split-exogenous")
        .arg("--input").arg(&bam)
        .arg("--output").arg(&prefix)
        .arg("--exogenous-prefix").arg("spike_")
        .arg("--allow-unknown-mapq");

    cmd.assert().success();

    for ext in ["endogenous.bam", "exogenous.bam", "both.bam", "unmapped.bam"] {
        let path = prefix.with_extension(ext);
        assert!(path.exists(), "missing output {ext}");
    }
}
