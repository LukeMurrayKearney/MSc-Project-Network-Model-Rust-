use crate::random_graphs::*;
// use crate::tau_leap::*;
use crate::run_model::*;
use crate::write_to_file::*;

pub fn test_run_model(n: usize) {
    // define network with initial infection
    let network_structure = NetworkStructure::new_molloy_reed(n, vec![n/9, 2*n/9, 3*n/9, 4*n/9, 5*n/9, 6*n/9, 7*n/9, 8*n/9, n], "");
    let mut network_properties = NetworkProperties::new(&network_structure);
    network_properties.params(vec![0.007, 3.0, 7.0, 1000.0]);
    network_properties.outbreak_type = OutbreakType::SEIRS;
    network_properties.result_type = ResultType::SecondaryCases(1);
    let maxtime: f64 = 365.0;
    let initially_infected = 0.001;

    let start = std::time::Instant::now();
    let output = run_model(&network_structure, &mut network_properties, maxtime, initially_infected);
    let elapsed = start.elapsed();
    println!("{} seconds", elapsed.as_secs());
    outbreak_results_csv(output, network_properties.result_type,"../../csv/test2.csv");
}

// pub fn test_sir_ba(n: usize) {

//     // define network with an initial infection
//     let initially_infected: f64 = 0.0001;
//     let network_structure: NetworkStructure = NetworkStructure::new_ba(n,5,3);
//     let mut network_properties: NetworkProperties = NetworkProperties::new(&network_structure);
//     network_properties.params(vec![0.02,0.2]);
    
//     // run SIR on network and time it 
//     let start = std::time::Instant::now();
//     let sir_results: Output = run_tau_leap(&network_structure, &mut network_properties, 5.0*365.0, 1.0, initially_infected);
//     let elapsed = start.elapsed();
//     println!("{} seconds", elapsed.as_secs());
//     outbreak_results_csv(sir_results, network_properties.result_type,"../../csv/BA_results_10_000.csv");
// }

// pub fn test_infecs_ba(n: usize, iters: usize) {

//     // define network and initial infection proportion
//     let initially_infected: f64 = 0.0001;
//     let network_structure: NetworkStructure = NetworkStructure::new_ba(n,5,3);
//     let mut network_properties: NetworkProperties = NetworkProperties::new(&network_structure);
    
//     // outbreak type and parameters
//     network_properties.params(vec![0.02, 1.0/7.0, 2.0/365.0]);
//     network_properties.result_type = ResultType::AvgInfections(iters);
    
//     // run outbreak with timing mechanism
//     let start = std::time::Instant::now();
//     let output = run_tau_parallel(&network_structure, &network_properties,5.0*365.0, 1.0, initially_infected);
//     let elapsed = start.elapsed();
//     println!("{} seconds", elapsed.as_secs());

//     //print results to a json
//     if let Err(err) = results_json(&output, "../../csv/BA_results_infections_SIRS_0.01.json") {
//         eprintln!("Error: {}", err);
//     }
// }

pub fn comix_sbm(n: usize) -> NetworkStructure {

    // age brackets, 0-17 / 18-39 / 40-65 / 65+
    // proportions,  20.7 / 28.5  / 32.2  / 19.6
    // https://www.ethnicity-facts-figures.service.gov.uk/uk-population-by-ethnicity/demographics/age-groups/latest
    let partitions: Vec<usize> = vec![
        (20.7*(n as f64)/100.0) as usize, 
        (49.2*(n as f64)/100.0) as usize, 
        (81.4*(n as f64)/100.0) as usize, 
        n
    ];

    let rates_mat: Vec<Vec<f64>> = vec![
        vec![5.3, 1.1, 2.1, 0.31],
        vec![1.1, 1.4, 2.2, 0.77],
        vec![2.1, 2.2, 2.4, 1.1],
        vec![0.31, 0.77, 1.1, 1.4]
    ];
    NetworkStructure::new_sbm(n, partitions, rates_mat)
}

pub fn comix_sbm_weighted(n: usize) -> NetworkStructure {

    // age brackets, 0-17 / 18-39 / 40-65 / 65+
    // proportions,  20.7 / 28.5  / 32.2  / 19.6nrows
    // https://www.ethnicity-facts-figures.service.gov.uk/uk-population-by-ethnicity/demographics/age-groups/latest
    let partitions: Vec<usize> = vec![
        (20.7*(n as f64)/100.0) as usize, 
        (49.2*(n as f64)/100.0) as usize, 
        (81.4*(n as f64)/100.0) as usize, 
        n
    ];

    let rates_mat: Vec<Vec<f64>> = vec![
        vec![5.3, 1.1, 2.1, 0.31],
        vec![1.1, 1.4, 2.2, 0.77],
        vec![2.1, 2.2, 2.4, 1.1],
        vec![0.31, 0.77, 1.1, 1.4]
    ];
    NetworkStructure::new_sbm_weighted(n, partitions, rates_mat)
}

pub fn network_structure_json(network_structure: &NetworkStructure) {

    let mut output: Output = Output::new();
    output.network_struct = SerializeableNetwork::from(network_structure);
    //print results to a json
    if let Err(err) = results_json(&output, "../../json/SBM_network_weight.json") {
        eprintln!("Error: {}", err);
    }
}