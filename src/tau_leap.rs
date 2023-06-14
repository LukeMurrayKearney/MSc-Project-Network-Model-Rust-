use crate::random_graphs::*;
extern crate nalgebra as na;
use na::{DMatrix, DVector};
 
pub fn run_tau_leap(network:&mut Network, parameters: &Vec<f64>, maxtime: f64, dt: f64) {
    for i in 0..((maxtime/dt) as usize) {
        step_tau_leap(network, parameters, dt);
    }
}

fn step_tau_leap(network:&mut Network, parameters: &Vec<f64>, dt: f64) {
    
}