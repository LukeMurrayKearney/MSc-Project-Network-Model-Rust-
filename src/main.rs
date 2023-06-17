use networks::random_graphs::*;
use networks::tau_leap::run_tau_leap;
use networks::write_to_file::outbreak_results_csv;
extern crate nalgebra as na;

pub(crate) fn main() {
    // define network with an initial infection
    let n: usize = 100_000;
    let initially_infected: f64 = 0.001;
    let mut network:Network = Network::new_ba(n,5,2);
    network.initialize_infection(initially_infected);
    
    // run SIR on network and time it 
    let start = std::time::Instant::now();
    run_tau_leap(&mut network, &vec![0.03,0.1], 5.0*365.0, 1.0);
    let elapsed = start.elapsed();
    println!("{} seconds", elapsed.as_secs());
    outbreak_results_csv(network);
}

