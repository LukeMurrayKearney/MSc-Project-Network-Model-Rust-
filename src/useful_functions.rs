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