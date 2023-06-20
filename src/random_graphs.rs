extern crate nalgebra as na;
use na::{DVector};
use nalgebra_sparse::{coo::CooMatrix, csr::CsrMatrix};
extern crate random_choice;
use self::random_choice::random_choice;
use rand::prelude::*;
use serde::Serialize;

#[derive(Clone,Debug)]
pub enum State {
    Susceptible,
    Infected,
    Recovered
}

#[derive(Clone)]
pub enum ResultType {
    SIR,
    AvgInfections(usize)
}

#[derive(Debug)]
pub struct NetworkStructure {
    pub adjacency_matrix: CsrMatrix<f64>,
    pub degree: DVector<f64>
}

#[derive(Clone)]
pub struct NetworkProperties {
    pub nodal_states: Vec<State>,
    pub results: Vec<Vec<usize>>,
    pub result_type: ResultType,
    pub parameters: Vec<f64>
}

#[derive(Debug,Serialize)]
pub struct Output {
    pub sir: Vec<Vec<usize>>,
    pub infections: Vec<Vec<usize>>
}


impl NetworkStructure {
    pub fn new_ba(n: usize, m0: usize, m: usize) -> NetworkStructure {
        //check dimensions correct
        if m0 < m {
            println!("m0 must be greater than m");
        }
        // start with complete sparse graph,
        let mut coo_mat:CooMatrix<f64> = CooMatrix::new(n,n);
        let mut degrees: Vec<f64> = vec![0.0;n];

        for i in 0..m0 {
            for j in 0..m0 {
                coo_mat.push(i,j, 1.0);
            }
        }
        for i in 0..m0 {
            degrees[i] += 1.0
        }
        let nodes: Vec<usize> = (0..n).collect();
        // implement BA algorithm:
        for i in m0..n {
            let choices = random_choice().random_choice_f64(&nodes, &degrees, m);
            for j in choices.into_iter() {
                coo_mat.push(i, *j, 1.0);
                coo_mat.push(*j, i, 1.0);
                degrees[i] += 1.0;
                degrees[*j] += 1.0;
            }
        }
        let matrix: CsrMatrix<f64> = CsrMatrix::from(&coo_mat);
        // result network struct with adjacency matrix
        NetworkStructure {
            adjacency_matrix: matrix,
            degree: DVector::from_vec(degrees),
        }
    }
}

impl NetworkProperties {
    pub fn new(network: &NetworkStructure) -> NetworkProperties {
        NetworkProperties { 
            nodal_states: vec![State::Susceptible; network.degree.len()],
            results: Vec::new(),
            result_type: ResultType::SIR,
            parameters: vec![0.1,0.2]
        }
    }

    pub fn params(&mut self, beta: f64, gamma: f64) {
        self.parameters = vec![beta,gamma];
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
        self.results.push(self.count_states());
    }

    pub fn count_states(&self) -> Vec<usize> {
        let mut result: Vec<usize> = vec![0; 3];
        for state in self.nodal_states.iter() {
            match state {
                State::Susceptible => result[0] += 1,
                State::Infected => result[1] += 1,
                State::Recovered => result[2] += 1
            }
        }
        result
    }
}

impl Output {
    pub fn new() -> Output {
        Output { sir: Vec::new(), infections: Vec::new() }
    }
}