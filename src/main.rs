use networks::random_graphs::*;
use networks::tau_leap::run_tau_leap;
use networks::write_to_file::outbreak_results_csv;
extern crate nalgebra as na;

pub(crate) fn main() {
    // define network with an initial infection
    let n: usize = 100_000;
    let initially_infected: f64 = 0.001;
    let network_structure: NetworkStructure = NetworkStructure::new_ba(n,5,2);
    let mut network_properties: NetworkProperties = NetworkProperties::new(&network_structure);
    network_properties.initialize_infection(initially_infected);
    network_properties.params(0.1,0.2);
    
    // run SIR on network and time it 
    let start = std::time::Instant::now();
    let sir_results: Output = run_tau_leap(&network_structure, &mut network_properties, 5.0*365.0, 1.0);
    let elapsed = start.elapsed();
    println!("{} seconds", elapsed.as_secs());
    outbreak_results_csv(sir_results, network_properties.result_type);
}

