use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct DistributionParameters {
    pub lambda: Vec<Vec<f64>>,
    pub p_geom: Vec<Vec<f64>>,
    pub p: Vec<Vec<f64>>
}

impl DistributionParameters {
    pub fn new() -> DistributionParameters {
        DistributionParameters { lambda: Vec::new(), p_geom: Vec::new(), p: Vec::new() }
    }
}

pub fn count_buckets(values: Vec<f64>) -> Vec<i32> {
    let mut buckets = vec![0; values.iter().map(|&x| x as usize).max().unwrap()];
    for i in values.iter() {
        buckets[*i as usize - 1] += 1; 
    }
    buckets
}

pub fn rates_to_probabilities(rates_mat: Vec<Vec<f64>>, partitions: &Vec<usize>) -> Vec<Vec<f64>> {
    
    // find consecutive group sizes to turn rates to probabilities
    let mut group_sizes: Vec<usize> = partitions
        .windows(2)
        .map(|pair| {
            pair[1] - pair[0]
        })
        .collect();
    group_sizes.insert(0,partitions[0]);
    
    // transform rates matrix to probability matrix 
    rates_mat
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter().map(|rate| {
                rate / (group_sizes[i] as f64)
            })
            .collect()
        })
        .collect()
}