extern crate nalgebra as na;
use na::{DMatrix, DVector};
extern crate random_choice;
use self::random_choice::random_choice;
use rand::prelude::*;

#[derive(Clone,Debug)]
pub enum State {
    Susceptible,
    Exposed,
    Infected,
    Recovered
}

#[derive(Debug)]
pub struct Network {
    pub adjacency_matrix: DMatrix<f64>,
    pub degree: DVector<f64>,
    pub nodal_states: Vec<State>,
    pub results: Vec<Vec<usize>>
}

impl Network {
    pub fn new_ba(n: usize, m0: usize, m: usize) -> Network {
        //check dimensions correct
        if m0 < m {
            println!("m0 must be greater than m");
        }
        // start with complete graph
        let mut matrix = DMatrix::from_fn(n, n,
            |r, c| if (m0 > r) & (m0 > c) { 1.0 } else {0.0});
        let mut degrees: Vec<f64> = vec![0.0;n];
        for i in 0..m0 {
            degrees[i] += 1.0
        }
        let nodes: Vec<usize> = (0..n).collect();
        // implement BA algorithm:
        for i in m0..n {
            let choices = random_choice().random_choice_f64(&nodes, &degrees, m);
            for j in choices.into_iter() {
                matrix[(i,*j)] += 1.0;
                matrix[(*j,i)] += 1.0;
                degrees[i] += 1.0;
                degrees[*j] += 1.0;
            }
        }
        // result network struct with adjacency matrix
        Network {
            adjacency_matrix: matrix,
            degree: DVector::from_vec(degrees),
            nodal_states: vec![State::Susceptible; n],
            results: vec![vec![n, 0, 0, 0]]
        }
    }

    pub fn initialize_infection(&mut self, proportion_of_population: f64) {
        let number_of_infecteds: usize = match proportion_of_population as usize {
            0..=1 => {
                ((self.nodal_states.len() as f64) * proportion_of_population) as usize
            },
            _ => {
                println!("The proportion infected must be between 0 and 1");
                0
            }
        };
        // define random number generator
        let mut rng = rand::thread_rng();
        // shuffle indices and choose
        let mut indices: Vec<usize> = (0..self.nodal_states.len()).collect();
        indices.shuffle(&mut rng);
        for i in 0..number_of_infecteds {
            self.nodal_states[indices[i]] = State::Infected
        }
    }

    pub fn count_states(&self) -> Vec<usize> {
        let mut result: Vec<usize> = vec![0; 4];
        for state in self.nodal_states.iter() {
            match state {
                State::Susceptible => result[0] += 1,
                State::Exposed => result[1] += 1,
                State::Infected => result[2] += 1,
                State::Recovered => result[3] += 1
            }
        }
        result
    }
}
