// use networks::random_graphs::*;
// use networks::tau_leap::*;
// use networks::write_to_file::*;
use networks::{run_scenarios::*, random_graphs::NetworkStructure};
extern crate nalgebra as na;

pub(crate) fn main() {
    
    let n: usize = 100;

    // function to test outbreak on BA graph and write SIR compartments to csv (SIR)
    // test_sir_ba(n);
    
    // function to test paralellisation of outbreaks and record infections in a csv (SIRS)
    // let iters: usize = 30;
    // test_infecs_ba(n, iters);

    // Testing SBM distribution
    // let network_structure: NetworkStructure = comix_sbm(n);
    // network_structure_json(&network_structure);

    // Testing SBM with weights
    let network_structure = comix_sbm_weighted(n);
}

