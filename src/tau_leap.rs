use crate::random_graphs::*;
extern crate nalgebra as na;
use na::{DVector};
use rand::{seq::SliceRandom, rngs::ThreadRng};
use rand_distr::{Poisson,Distribution};
extern crate random_choice;
use self::random_choice::random_choice;

 
pub fn run_tau_leap(network:&mut Network, parameters: &Vec<f64>, maxtime: f64, dt: f64) {
    let mut rng = rand::thread_rng();
    for _ in 0..((maxtime/dt) as usize) {
        step_tau_leap(network, parameters, dt, &mut rng);
    }
}

fn step_tau_leap(network:&mut Network, parameters: &Vec<f64>, dt: f64, rng: &mut ThreadRng) {
    // define vector of infected and susceptible individuals
    let infecteds: DVector<f64> = DVector::from_fn(network.nodal_states.len(), 
        |i, _| match network.nodal_states[i] {
            State::Infected => 1.0,
            _ => 0.0,
        }
    );
    let susceptibles: DVector<f64> = DVector::from_fn(network.nodal_states.len(), 
        |i, _| match network.nodal_states[i] {
            State::Susceptible => 1.0,
            _ => 0.0,
        }
    );
    // calculate rates and poisson samples
    let rates: (Vec<f64>, f64) = calculate_rates(network, infecteds.clone(), parameters);
    let total_rate: Vec<f64> = vec![rates.0.iter().sum(), rates.1];
    let num_events: Vec<usize> = poisson_rvs(total_rate, dt, rng);
    update_compartments(network, num_events, susceptibles, infecteds, rates.0, rng);
    network.results.push(network.count_states());
}

fn calculate_rates(network:&mut Network, infecteds: DVector<f64>, parameters: &Vec<f64>) -> (Vec<f64>, f64) {
    // rates of each occurence 
    let ds_vec = parameters[0] * (infecteds.transpose() * network.adjacency_matrix.clone());
    let di = parameters[1] * infecteds.sum();
    let ds: Vec<f64> = ds_vec
        .iter()
        .enumerate()
        .filter(|&(index, _)| 
            match network.nodal_states[index] {
                State::Susceptible => true,
                _ => false
            }
        )
        .map(|(_, x)| *x)
        .collect();
    (ds,di)
}

fn poisson_rvs(rates: Vec<f64>, dt: f64, rng: &mut ThreadRng) -> Vec<usize> {
    // create a random number generator, and generate poisson distributed random numbers with a given rate
    let mut result: Vec<usize> = Vec::new(); 
    for lambda in rates.iter() {
        let poi = Poisson::new(*lambda * dt).unwrap();
        result.push(poi.sample(rng) as usize);
    }
    result
}

fn update_compartments(network: &mut Network, num: Vec<usize>, susceptibles: DVector<f64>, infecteds: DVector<f64>, 
    infection_pressure: Vec<f64>, rng: &mut ThreadRng) {

    let mut x: usize;
    // check values are possible 
    if num[0] < network.results.last().unwrap()[0] {
        x = num[0];
    } else {
        x = network.results.last().unwrap()[0];
    }
    // get susceptible indices
    let indices: Vec<usize> = susceptibles
        .iter()
        .enumerate()
        .filter(|(_, &value)| value != 0.0)
        .map(|(index, _)| index)
        .collect();
    // chooose people to infect and update
    let choices: Vec<&usize> = random_choice().random_choice_f64(&indices, &infection_pressure, x);
    for i in choices.iter() {
        network.nodal_states[**i] = State::Infected;
    }

    // Repeat for recovery events
    if num[1] < network.results.last().unwrap()[2] {
        x = num[1];
    } else {
        x = network.results.last().unwrap()[2];
    }
    let mut indices: Vec<usize> = infecteds
        .iter()
        .enumerate()
        .filter(|(_, &value)| value != 0.0)
        .map(|(index, _)| index)
        .collect();
    indices.shuffle(rng);
    for i in 0..x {
        network.nodal_states[indices[i]] = State::Recovered;
    }
}