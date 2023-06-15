use networks::random_graphs::*;
use networks::tau_leap::*;
extern crate nalgebra as na;
use na::{DMatrix, DVector};

pub(crate) fn main() {
    // define network with an initial infection
    let n: usize = 100;
    let initially_infected: f64 = 0.1;
    let mut network:Network = Network::new_ba(n,3usize,3usize);
    network.initialize_infection(initially_infected);
    
    // run SIR on network
    
    println!("{:?}, \n {:?}", network.nodal_states, network.count_states());
    run_tau_leap(&mut network, &vec![0.1,0.1], 1.0, 1.0);
}

