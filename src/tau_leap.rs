use std::iter::Sum;

use crate::random_graphs::*;
extern crate nalgebra as na;
use na::{DMatrix, DVector};
use ndarray::iter::IterMut;
use rand_distr::{Poisson,Distribution};

 
pub fn run_tau_leap(network:&mut Network, parameters: &Vec<f64>, maxtime: f64, dt: f64) {
    for i in 0..((maxtime/dt) as usize) {
        step_tau_leap(network, parameters, dt);
    }
}

fn step_tau_leap(network:&mut Network, parameters: &Vec<f64>, dt: f64) {
    // define vector of infected individuals
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
    let rates: Vec<f64> = calculate_rates(network, infecteds, parameters);
    let mut num_instances: Vec<usize> = poisson_rvs(rates);
    
    
}

fn calculate_rates(network:&mut Network, infecteds: DVector<f64>, parameters: &Vec<f64>) -> Vec<f64> {
    // rates of each occurence 
    let ds_vec = parameters[0] * (infecteds.transpose() * network.adjacency_matrix.clone());
    let di = parameters[1] * infecteds.sum();
    let ds: f64 = ds_vec
        .iter()
        .enumerate()
        .filter(|&(index, _)| 
            match network.nodal_states[index] {
                State::Susceptible => true,
                _ => false
            }
        )
        .map(|(_, x)| *x)
        .sum();
    vec![ds,di]
}

fn poisson_rvs(rates: Vec<f64>) -> Vec<usize> {
    // create a random number generator, and generate poisson distributed random numbers with a given rate
    let mut rng = rand::thread_rng();
    let mut result: Vec<usize> = Vec::new(); 
    for lambda in rates.iter() {
        let poi = Poisson::new(*lambda).unwrap();
        result.push(poi.sample(&mut rng) as usize);
    }
    result
}

fn update_compartments(network: &mut Network, num: Vec<usize>, susceptibles: DVector<f64>, infecteds: DVector<f64>) {
    let x: usize;
    // check values are possible and update compartments
    if num[0] < network.results.last().unwrap()[0] {
        x = num[0]
    } else {
        x = network.results.last().unwrap()[0];
    }
    
}