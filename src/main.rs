use networks::random_graphs::*;
use networks::tau_leap::*;
extern crate nalgebra as na;
// use na::{DMatrix, DVector};

pub(crate) fn main() {
    // define network with an initial infection
    let n: usize = 10_000;
    let initially_infected: f64 = 0.001;
    let mut network:Network = Network::new_ba(n,10usize,2usize);
    network.initialize_infection(initially_infected);
    
    // run SIR on network
    run_tau_leap(&mut network, &vec![0.05,0.1], 100.0, 1.0);
}

