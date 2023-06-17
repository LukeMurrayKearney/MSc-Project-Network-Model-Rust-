use crate::random_graphs::*;
extern crate nalgebra as na;
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use rand::{seq::{SliceRandom}, rngs::ThreadRng};
use rand_distr::{Poisson,Distribution};
extern crate random_choice;
use self::random_choice::random_choice;

 
pub fn run_tau_leap(network:&mut Network, parameters: &Vec<f64>, maxtime: f64, dt: f64) {
    let mut rng = rand::thread_rng();
    for i in 0..((maxtime/dt) as usize) {
        step_tau_leap(network, parameters, dt, &mut rng);
        if i % 10 == 0 {
            println!("{i}");
        }
        if network.results.last().unwrap()[1] == 0 {
            break;
        }
    }
}

fn step_tau_leap(network:&mut Network, parameters: &Vec<f64>, dt: f64, rng: &mut ThreadRng) {
    let mut sir_coo: (CooMatrix<f64>, CooMatrix<f64>, CooMatrix<f64>) = 
        (
            CooMatrix::new(1, network.nodal_states.len()),
            CooMatrix::new(1, network.nodal_states.len()),
            CooMatrix::new(1, network.nodal_states.len()),
        );

    for (i, x) in network.nodal_states.iter().enumerate() {
        match x {
            State::Susceptible => sir_coo.0.push(0, i, 1.0),
            State::Infected => sir_coo.1.push(0, i, 1.0),
            State::Recovered => sir_coo.2.push(0, i, 1.0)
        }
    } 
    // calculate rates and poisson samples
    let rates: (Vec<f64>, f64) = calculate_rates(network, parameters, &sir_coo.1, &sir_coo.0);
    let total_rate: Vec<f64> = vec![rates.0.iter().sum(), rates.1];
    let num_events: Vec<usize> = poisson_rvs(total_rate, dt, rng);
    update_compartments(network, num_events, sir_coo.0, sir_coo.1, rates.0, rng);
    network.results.push(network.count_states());
}

fn calculate_rates(network:&mut Network, parameters: &Vec<f64>, infecteds_coo: &CooMatrix<f64>, susceptibles_coo: &CooMatrix<f64>) -> (Vec<f64>, f64) {
    // rates of each occurence 
    let infecteds: CsrMatrix<f64> = CsrMatrix::from(infecteds_coo);
    let ds_vec = parameters[0] * (infecteds.clone() * network.adjacency_matrix.clone());
    let total_infecteds: f64 = infecteds.values().iter().sum();
    let di = parameters[1] * total_infecteds;
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

fn update_compartments(network: &mut Network, num: Vec<usize>, susceptibles: CooMatrix<f64>, infecteds: CooMatrix<f64>, 
    infection_pressure: Vec<f64>, rng: &mut ThreadRng) {

    let mut x: usize;
    // check values are possible 
    if num[0] < network.results.last().unwrap()[0] {
        x = num[0];
    } else {
        x = network.results.last().unwrap()[0];
    }
    // get susceptible indices
    let indices: Vec<usize> = susceptibles.col_indices().to_owned();
    // chooose people to infect and update
    let choices: Vec<&usize> = random_choice().random_choice_f64(&indices, &infection_pressure, x);
    for i in choices.iter() {
        network.nodal_states[**i] = State::Infected;
    }

    // Repeat for recovery events
    if num[1] < network.results.last().unwrap()[1] {
        x = num[1];
    } else {
        x = network.results.last().unwrap()[1];
    }
    let mut indices: Vec<usize> = infecteds.col_indices().to_owned();
    indices.shuffle(rng);
    for i in 0..x {
        network.nodal_states[indices[i]] = State::Recovered;
    }
}