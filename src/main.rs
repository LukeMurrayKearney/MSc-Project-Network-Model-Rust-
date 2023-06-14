use networks::random_graphs::*;
extern crate nalgebra as na;
use na::{DMatrix, DVector};

pub(crate) fn main() {
    // define network
    let n: usize = 10;
    let mut network:Network = Network::new_ba(n,3usize,3usize);

    // run SIR on network
    network.nodal_states = DVector::zeros(n);
    println!("{:?}, \n {:?}", network.degree, network.nodal_states);
}

