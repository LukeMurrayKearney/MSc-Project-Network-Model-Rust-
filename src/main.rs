use networks::random_graphs::*;
use networks::run_model::*;
// use networks::tau_leap::*;
use networks::write_to_file::*;
use networks::{run_scenarios::*, random_graphs::NetworkStructure, useful_functions::*};
extern crate nalgebra as na;

pub(crate) fn main() {
    
    // Network size
    let n: usize = 20_000;

    // Testing SBM
    // let network_structure: NetworkStructure = comix_sbm_weighted(n);
    // network_structure_json(&network_structure);

    // Outputting degree distribution for molloy reed
    // let network = NetworkStructure::new_molloy_reed(n, vec![n/9, 2*n/9, 3*n/9, 4*n/9, 5*n/9, 6*n/9, 7*n/9, 8*n/9, n], "");
    // let count = count_buckets(network.degree)
    //     .iter()
    //     .enumerate()
    //     .map(|(i,x)| (i, *x as usize))
    //     .collect();
    // let result = vector_to_csv(count, "model_output_files/MR_dd_100k.csv");

    // Testing new model 
    test_run_model(n);
}

