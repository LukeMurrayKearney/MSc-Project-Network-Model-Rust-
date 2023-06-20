// use std::thread::Result;
use std::sync::{Arc, Mutex};
use crate::random_graphs::*;
extern crate nalgebra as na;
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use rand::{seq::{SliceRandom}, rngs::ThreadRng};
use rand_distr::{Poisson,Distribution};
extern crate random_choice;
use self::random_choice::random_choice;
use rayon::prelude::*;

pub fn run_tau_parallel(network_structure: &NetworkStructure, network_properties: &NetworkProperties, maxtime: f64, dt: f64, initially_infected: f64) -> Output {
    let iterations: Vec<usize> = match network_properties.result_type {
        ResultType::SIR => vec![0],
        ResultType::AvgInfections(iters) => (0..iters).collect()
    };
    let output = Arc::new(Mutex::new(Output::new()));
    iterations.par_iter().for_each(|i| {
        let tmp = run_tau_leap(network_structure, &mut network_properties.clone(), maxtime, dt, initially_infected);
        match network_properties.result_type {
            ResultType::SIR => {
                let mut output = output.lock().unwrap();
                output.sir = tmp.sir;
            },
            ResultType::AvgInfections(_) => {
                let mut output = output.lock().unwrap();
                output.infections.push(tmp.infections[0].clone());
            }
        }
        println!("{i}");
    });
    Arc::try_unwrap(output).unwrap().into_inner().unwrap() // Unwrap the `Mutex` and return the `Output`
}

pub fn run_tau_leap(network_structure: &NetworkStructure, network_properties: &mut NetworkProperties, maxtime: f64, dt: f64, initially_infected: f64) -> Output {
    network_properties.initialize_infection(initially_infected);
    let mut rng = rand::thread_rng();
    for _i in 0..((maxtime/dt) as usize) {
        step_tau_leap(&network_structure, network_properties, dt, &mut rng);
        // if i % 100 == 0 {
        //     println!("{i}");
        // }
        if network_properties.results.last().unwrap()[1] == 0 {
            break;
        }
    }
    // matching the measures wanted from
    let mut output: Output = Output::new(); 
    match network_properties.result_type {
        ResultType::SIR => {
            output.sir = network_properties.results.clone();
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
        }
    }
    output
}

fn step_tau_leap(network_structure: &NetworkStructure, network_properties: &mut NetworkProperties, dt: f64, rng: &mut ThreadRng) {
    let mut sir_coo: (CooMatrix<f64>, CooMatrix<f64>, CooMatrix<f64>) = 
        (
            CooMatrix::new(1, network_properties.nodal_states.len()),
            CooMatrix::new(1, network_properties.nodal_states.len()),
            CooMatrix::new(1, network_properties.nodal_states.len()),
        );

    for (i, x) in network_properties.nodal_states.iter().enumerate() {
        match x {
            State::Susceptible => sir_coo.0.push(0, i, 1.0),
            State::Infected => sir_coo.1.push(0, i, 1.0),
            State::Recovered => sir_coo.2.push(0, i, 1.0)
        }
    } 
    // calculate rates and poisson samples
    let rates: (Vec<f64>, f64) = calculate_rates(network_structure, network_properties, &sir_coo.1, &sir_coo.0);
    let total_rate: Vec<f64> = vec![rates.0.iter().sum(), rates.1];
    let num_events: Vec<usize> = poisson_rvs(total_rate, dt, rng);
    update_compartments(network_properties, num_events, sir_coo.0, sir_coo.1, rates.0, rng);
    network_properties.results.push(network_properties.count_states());
}

fn calculate_rates(network_structure: &NetworkStructure, network_properties: &mut NetworkProperties, infecteds_coo: &CooMatrix<f64>, susceptibles_coo: &CooMatrix<f64>) -> (Vec<f64>, f64) {
    // rates of each occurence 
    let infecteds: CsrMatrix<f64> = CsrMatrix::from(infecteds_coo);
    let ds_vec = network_properties.parameters[0] * (infecteds.clone() * network_structure.adjacency_matrix.clone());
    let total_infecteds: f64 = infecteds.values().iter().sum();
    let di = network_properties.parameters[1] * total_infecteds;
    //must find common indices to apply rates to the susceptible individuals only
    let indices_infection_rate: Vec<usize> = ds_vec.col_indices().to_owned();
    let indices_susceptibles: Vec<usize> = susceptibles_coo.col_indices().to_owned();
    let common_indices: Vec<usize> = indices_infection_rate
        .iter()
        .enumerate()
        .filter(|(_, &x)| indices_susceptibles.contains(&x))
        .map(|(i, _)| i)
        .collect();
    // with common indices reduce infection rate vector to only susceptible individuals
    let mut ds: Vec<f64> = Vec::new();
    for index in common_indices.iter() {
        ds.push(ds_vec.values()[*index]); 
        
    }
    (ds,di)
}

fn poisson_rvs(rates: Vec<f64>, dt: f64, rng: &mut ThreadRng) -> Vec<usize> {
    // create a random number generator, and generate poisson distributed random numbers with a given rate
    let mut result: Vec<usize> = Vec::new(); 
    for lambda in rates.iter() {
        if *lambda > std::f64::EPSILON {
            let poi = Poisson::new(*lambda * dt).unwrap();
            result.push(poi.sample(rng) as usize);
        } else {
            result.push(0)
        }
    }
    result
}

fn update_compartments(network_properties: &mut NetworkProperties, num: Vec<usize>, susceptibles: CooMatrix<f64>, infecteds: CooMatrix<f64>, 
    infection_pressure: Vec<f64>, rng: &mut ThreadRng) {

    let mut x: usize;
    // check values are possible 
    if num[0] < network_properties.results.last().unwrap()[0] {
        x = num[0];
    } else {
        x = network_properties.results.last().unwrap()[0];
    }
    // get susceptible indices
    let indices: Vec<usize> = susceptibles.col_indices().to_owned();
    // chooose people to infect and update
    let choices: Vec<&usize> = random_choice().random_choice_f64(&indices, &infection_pressure, x);
    for i in choices.iter() {
        network_properties.nodal_states[**i] = State::Infected;
    }

    // Repeat for recovery events
    if num[1] < network_properties.results.last().unwrap()[1] + x {
        x = num[1];
    } else {
        x += network_properties.results.last().unwrap()[1];
    }
    let mut indices: Vec<usize> = infecteds.col_indices().to_owned();
    // add new infection events to selection
    for i in choices.into_iter() {
        indices.push(*i)
    }
    indices.shuffle(rng);
    //choose people to recover
    for i in 0..x {
        network_properties.nodal_states[indices[i]] = State::Recovered;
    }
}