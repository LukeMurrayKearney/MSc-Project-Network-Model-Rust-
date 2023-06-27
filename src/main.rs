use std::vec;

use na::DMatrix;
// use networks::random_graphs::*;
// use networks::tau_leap::*;
// use networks::write_to_file::*;
use networks::{run_scenarios::*, random_graphs::NetworkStructure};
extern crate nalgebra as na;

pub(crate) fn main() {
    
    let n: usize = 4;

    // function to test outbreak on BA graph and write SIR compartments to csv (SIR)
    // test_sir_ba(n);
    
    // function to test paralellisation of outbreaks and record infections in a csv (SIRS)
    // let iters: usize = 30;
    // test_infecs_ba(n, iters);

    // Testing SBM
    let partitions: Vec<usize> = vec![n/2,n];
    let rates_mat: Vec<Vec<f64>> = vec![
        vec![1.0,0.0],
        vec![0.0,1.0]
    ];
    let network_structure: NetworkStructure = NetworkStructure::new_sbm(n, partitions, rates_mat);
    let matrix = DMatrix::from(&network_structure.adjacency_matrix);
    dbg!(matrix);
}

