use anyhow::Result;
use polars::prelude::*;

/// Collapse *adjacent* bins (contiguous end==next start) with identical `score` into single intervals.
pub fn collapse_adjacent_bins(df: DataFrame) -> Result<DataFrame> {
    if df.height() == 0 {
        return Ok(df);
    }

    let lf = df
        .lazy()
        .sort(["chrom", "start"], Default::default())
        .with_columns([
            col("chrom").shift(lit(1)).alias("prev_chrom"),
            col("end").shift(lit(1)).alias("prev_end"),
            col("score").shift(lit(1)).alias("prev_score"),
        ]);

    // Start a new group if:
    // - first row, or
    // - chromosome changes, or
    // - bins are not contiguous (prev_end != start), or
    // - score changes (within epsilon)
    let eps = 1e-6_f64;
    let score_changed = (col("score") - col("prev_score")).abs().gt(lit(eps));

    let new_group = col("prev_chrom")
        .is_null()
        .or(col("chrom").neq(col("prev_chrom")))
        .or(col("start").neq(col("prev_end")))
        .or(score_changed);

    let collapsed = lf
        .with_column(new_group.cast(DataType::UInt32).cum_sum(false).alias("grp"))
        .group_by_stable([col("chrom"), col("grp")])
        .agg([
            col("start").min().alias("start"),
            col("end").max().alias("end"),
            col("score").first().alias("score"),
        ])
        .select([col("chrom"), col("start"), col("end"), col("score")])
        .collect()?;

    Ok(collapsed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn df_rows(df: &DataFrame) -> Result<Vec<(String, i64, i64, f64)>> {
        let chrom = df.column("chrom")?.str()?;
        let start = df.column("start")?.i64()?;
        let end = df.column("end")?.i64()?;
        let score = df.column("score")?.f64()?;
        let mut rows = Vec::with_capacity(df.height());
        for i in 0..df.height() {
            rows.push((
                chrom.get(i).unwrap().to_string(),
                start.get(i).unwrap(),
                end.get(i).unwrap(),
                score.get(i).unwrap(),
            ));
        }
        Ok(rows)
    }

    #[test]
    fn test_collapse_empty_df() -> Result<()> {
        let df = DataFrame::new(vec![
            Column::new("chrom".into(), Vec::<&str>::new()),
            Column::new("start".into(), Vec::<i64>::new()),
            Column::new("end".into(), Vec::<i64>::new()),
            Column::new("score".into(), Vec::<f64>::new()),
        ])?;
        let out = collapse_adjacent_bins(df)?;
        assert_eq!(out.height(), 0);
        Ok(())
    }

    #[test]
    fn test_collapse_contiguity_and_chrom_scoping() -> Result<()> {
        let df = DataFrame::new(vec![
            Column::new(
                "chrom".into(),
                &["chr1", "chr1", "chr1", "chr1", "chr2", "chr2"],
            ),
            Column::new("start".into(), &[0i64, 10, 21, 30, 0, 5]),
            Column::new("end".into(), &[10i64, 20, 30, 40, 5, 10]),
            Column::new("score".into(), &[1.0f64, 1.0, 1.0, 2.0, 1.0, 1.0]),
        ])?;

        let collapsed = collapse_adjacent_bins(df)?;
        let rows = df_rows(&collapsed)?;

        let expected = vec![
            ("chr1".to_string(), 0, 20, 1.0),
            ("chr1".to_string(), 21, 30, 1.0),
            ("chr1".to_string(), 30, 40, 2.0),
            ("chr2".to_string(), 0, 10, 1.0),
        ];

        assert_eq!(rows, expected);
        Ok(())
    }

    #[test]
    fn test_collapse_within_epsilon() -> Result<()> {
        let df = DataFrame::new(vec![
            Column::new("chrom".into(), &["chr1", "chr1", "chr1"]),
            Column::new("start".into(), &[0i64, 10, 20]),
            Column::new("end".into(), &[10i64, 20, 30]),
            Column::new("score".into(), &[1.0f64, 1.0000004, 1.000002]),
        ])?;

        let collapsed = collapse_adjacent_bins(df)?;
        let rows = df_rows(&collapsed)?;

        let expected = vec![
            ("chr1".to_string(), 0, 20, 1.0),
            ("chr1".to_string(), 20, 30, 1.000002),
        ];

        assert_eq!(rows, expected);
        Ok(())
    }
}
