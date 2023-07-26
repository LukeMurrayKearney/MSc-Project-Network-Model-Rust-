use crate::random_graphs::*;
use rand::{rngs::ThreadRng, Rng};
use rand_distr::{Distribution, Poisson};

// pub fn run_model_parallel(network_structure: &NetworkStructure, network_properties: &NetworkProperties, maxtime: f64, dt: f64, initially_infected: f64) -> Output {
//     let iterations: Vec<usize> = match network_properties.result_type {
//         ResultType::SIR => vec![0],
//         ResultType::AvgInfections(iters) => (0..iters).collect()
//     };
//     let output = Arc::new(Mutex::new(Output::new()));
//     iterations.par_iter().for_each(|i| {
//         let tmp = run_tau_leap(network_structure, &mut network_properties.clone(), maxtime, dt, initially_infected);
//         match network_properties.result_type {
//             ResultType::SIR => {
//                 let mut output = output.lock().unwrap();
//                 output.sir = tmp.sir;
//             },
//             ResultType::AvgInfections(_) => {
//                 let mut output = output.lock().unwrap();
//                 output.infections.push(tmp.infections[0].clone());
//             }
//         }
//         println!("{i}");
//     });
//     Arc::try_unwrap(output).unwrap().into_inner().unwrap() // Unwrap the `Mutex` and return the `Output`
// }

pub fn run_model(network_structure: &NetworkStructure, network_properties: &mut NetworkProperties, maxtime: f64, initially_infected: f64) -> Output {
    network_properties.initialize_infection(initially_infected);
    let mut rng = rand::thread_rng();
    for i in 0..(maxtime as usize) {
        step_model(&network_structure, network_properties, &mut rng);
        if i % 10 == 0 {
            println!("{i}");
        }
        if network_properties.results.last().unwrap()[1] + network_properties.results.last().unwrap()[2] == 0 {
            break;
        }
    }
    // matching the measures wanted from
    let mut output: Output = Output::new(); 
    match network_properties.result_type {
        ResultType::SEIR => {
            output.seir = network_properties.results.clone();
        },
        ResultType::AvgInfections(_) => {
            let mut infections: Vec<Vec<usize>> = Vec::new();
            infections.push(
                network_properties.results
                    .iter()
                    .map(|x| x[1])
                    .collect()
            );
            output.infections = infections;
        },
        ResultType::SecondaryCases(_) => {
            output.secondary_cases.push(
                network_properties.secondary_cases
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| {
                        match network_properties.nodal_states[*i] {
                            State::Recovered(_) => true,
                            _ => false
                        }
                    })
                    .map(|(_, x)| *x)
                    .collect()
            );
        }
    }
    output
}

fn step_model(network_structure: &NetworkStructure, network_properties: &mut NetworkProperties, rng: &mut ThreadRng) {
    let mut next_states: Vec<State> = vec![State::Susceptible; network_structure.degree.len()];
    let poisson_infectious = Poisson::new(network_properties.parameters[2]).unwrap();
    for (i, state) in network_properties.nodal_states.iter().enumerate() {
        match *state {
            State::Susceptible => (),
            State::Exposed(days) => {
                if days > network_properties.parameters[1] as usize {
                    next_states[i] = State::Infected(poisson_infectious.sample(rng) as usize);
                }
                else {
                    next_states[i] = State::Exposed(days + 1);
                }
            },
            State::Infected(days) => {
                if days <= 0 {
                    next_states[i] = State::Recovered(0);
                }
                else {
                    next_states[i] = State::Infected(days - 1);
                }
                // find connections to infected individuals
                let connections: Vec<usize> = network_structure.adjacency_matrix
                    .triplet_iter()
                    .filter(|&(row, _, _)| row == i)
                    .map(|(_, j, _)| j)
                    .collect();
                for j in connections.iter() {
                    match network_properties.nodal_states[*j] {
                        State::Susceptible => {
                            if rng.gen::<f64>() < network_properties.parameters[0] {
                                next_states[*j] = State::Exposed(0);
                                network_properties.secondary_cases[i] += 1;
                            }
                        },
                        _ => ()
                    }
                }
            },
            State::Recovered(days) => {
                if days < network_properties.parameters[3] as usize {
                    next_states[i] = State::Recovered(days + 1);
                }
            }
        }
    }
    network_properties.nodal_states = next_states;
    network_properties.results.push(network_properties.count_states());
}
