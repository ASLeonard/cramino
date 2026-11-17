/// Computes the N-statistic (e.g. N50 when `percentile` is 0.50) from lengths sorted descending.
pub fn get_n(lengths: &[u32], nb_bases_total: u64, percentile: f64) -> u64 {
    if lengths.is_empty() {
        return 0;
    }
    let mut acc: u64 = 0;
    for &val in lengths {
        acc += u64::from(val);
        if acc as f64 > nb_bases_total as f64 * percentile {
            return u64::from(val);
        }
    }

    u64::from(lengths[lengths.len() - 1])
}

/// Computes the Phred-scaled mean accuracy (QV) from the sum of per-read error rates.
pub fn mean_accuracy(qual_sum: f64, qual_reads: usize) -> f64 {
    if qual_reads == 0 {
        return f64::NAN;
    }
    let mean_err = qual_sum / qual_reads as f64;
    if mean_err <= 1e-6 {
        60.0
    } else {
        -10.0 * mean_err.log10()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_n_and_mean_accuracy_handle_normal_and_edge_cases() {
        assert_eq!(get_n(&[], 0, 0.50), 0);
        assert!(mean_accuracy(0.0, 0).is_nan());
        assert_eq!(mean_accuracy(0.0, 5), 60.0);
        assert_eq!(mean_accuracy(1e-7, 1), 60.0);

        let lengths = vec![100, 50, 25, 25];
        assert_eq!(get_n(&lengths, 200, 0.50), 50);

        assert!((mean_accuracy(0.02, 2) - 20.0).abs() < 1e-9);
    }
}
