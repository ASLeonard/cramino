
pub fn get_n(lengths: &Vec<u128>, nb_bases_total: u128, percentile: f64) -> u128 {
    let mut acc = 0;
    for val in lengths.iter() {
        acc += *val;
        if acc as f64 > nb_bases_total as f64 * percentile {
            return *val;
        }
    }

    lengths[lengths.len() - 1]
}

pub fn mean_accuracy(qualities: &Vec<f64>) -> f64 {
    -10.0 * (qualities.iter().sum::<f64>() / qualities.len() as f64).log10()
}
