use networks::random_graphs::*;
// use networks::tau_leap::*;
use networks::write_to_file::*;
use networks::{run_scenarios::*, random_graphs::NetworkStructure};
extern crate nalgebra as na;

pub(crate) fn main() {
    
    let n: usize = 10;

    // function to test outbreak on BA graph and write SIR compartments to csv (SIR)
    // test_sir_ba(n);
    
    // function to test paralellisation of outbreaks and record infections in a csv (SIRS)
    // let iters: usize = 30;
    // test_infecs_ba(n, iters);

    // Testing SBM
    // let network_structure: NetworkStructure = comix_sbm_weighted(n);
    // network_structure_json(&network_structure);

    // reading in rates matrix
    let network = NetworkStructure::new_molloy_reed(n, vec![5,n]);
    dbg!(network);
}

