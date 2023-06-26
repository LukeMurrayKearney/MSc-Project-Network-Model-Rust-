// use networks::random_graphs::*;
// use networks::tau_leap::*;
// use networks::write_to_file::*;
use networks::run_scenarios::*;
extern crate nalgebra as na;


pub(crate) fn main() {
    
    let n: usize = 100_000;

    // function to test outbreak on BA graph and write SIR compartments to csv (SIR)
    test_sir_ba(n);
    
    // function to test paralellisation of outbreaks and record infections in a csv (SIRS)
    test_infecs_ba(n);

    // Testing SBM

}

