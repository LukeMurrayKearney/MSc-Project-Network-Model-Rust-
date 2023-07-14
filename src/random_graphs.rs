use crate::useful_functions::*;
use crate::write_to_file::read_params_json;
extern crate nalgebra as na;
use std::vec;
use nalgebra_sparse::{coo::CooMatrix, csr::CsrMatrix};
extern crate random_choice;
use self::random_choice::random_choice;
use statrs::distribution::{Poisson, Geometric, NegativeBinomial};
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

#[derive(Clone)]
pub enum OutbreakType {
    SIR,
    SIRS
}

#[derive(Debug)]
pub struct NetworkStructure {
    pub adjacency_matrix: CsrMatrix<f64>,
    pub degree: Vec<f64>,
    pub age_brackets: Vec<usize>
}

#[derive(Clone)]
pub struct NetworkProperties {
    pub nodal_states: Vec<State>,
    pub results: Vec<Vec<usize>>,
    pub result_type: ResultType,
    pub outbreak_type: OutbreakType,
    pub parameters: Vec<f64>
}

#[derive(Debug,Serialize)]
pub struct Output {
    pub sir: Vec<Vec<usize>>,
    pub infections: Vec<Vec<usize>>,
    pub network_struct: SerializeableNetwork
}

impl NetworkStructure {

    pub fn new_molloy_reed(n: usize, partitions: Vec<usize>, file_path: &str) -> NetworkStructure {     
        let mut rng: ThreadRng = rand::thread_rng();
        let mut coo_mat: CooMatrix<f64> = CooMatrix::new(n,n);
        let mut degrees: Vec<f64> = vec![0.0;n];
        let mut count: usize = 0;
        // calculate group sizes
        let mut group_sizes: Vec<usize> = partitions
            .windows(2)
            .map(|pair| {
                pair[1] - pair[0]
            })
            .collect();
        group_sizes.insert(0,partitions[0]);
        // import parameters to sample
        let dist_params = read_params_json(file_path);
        
        // start iteration through age brackets, take(i+1) makes loop run over lower diag to remove double counting
        for (i, x) in partitions.iter().enumerate() {
            for (j, y) in partitions.iter().enumerate().take(i+1) {
                // create distributions and sample degrees from this 
                let poisson = Poisson::new(dist_params.lambda[i][j]);
                let geometric = Geometric::new(dist_params.p_geom[i][j]);
                let p = dist_params.p[i][j];
                let mut out_degree: Vec<usize> = (0..group_sizes[i])
                    .map(|_| {
                        (p*poisson.as_ref().unwrap().sample(&mut rng) + (1.0-p)*geometric.as_ref().unwrap().sample(&mut rng)) as usize
                    }
                    )
                    .collect();
                // repeat for in degrees
                let poisson = Poisson::new(dist_params.lambda[j][i]);
                let geometric = Geometric::new(dist_params.p_geom[j][i]);
                let p = dist_params.p[j][i];
                let mut in_degree: Vec<usize> = (0..group_sizes[j])
                    .map(|_| {
                        (p*poisson.as_ref().unwrap().sample(&mut rng) + (1.0-p)*geometric.as_ref().unwrap().sample(&mut rng)) as usize
                    }
                    )
                    .collect();
                
                // loop over out degrees and match with the in nodes
                let mut repeats: Vec<(usize,usize)> = Vec::new();
                for node_i in 0..out_degree.len() {
                    // don't allow self loops and double counting
                    if i == j {
                        in_degree[node_i] = 0;
                    }
                    count = 0;
                    let mut connections: Vec<usize> = in_degree
                        .iter()
                        .enumerate()
                        .filter(|(idx, x)| {
                            // stop if passed out degree, and in degree must be non-zero
                            if count < out_degree[node_i] && **x > 0 {
                                // check if link already created
                                if *idx >= node_i {
                                    count += 1;
                                    return true
                                }
                                else if !repeats.contains(&(*idx, node_i)) {
                                    count += 1;
                                    return true
                                }
                                else {
                                    return false
                                }
                            }
                            return false
                        })
                        .map(|(i,_)| i)
                        .collect();
                    connections.shuffle(&mut rng);
                    // connect nodes in the adjacency matrix
                    for node_j in connections.iter() {
                        coo_mat.push(*x-group_sizes[i]+node_i, *y-group_sizes[j]+*node_j, 1.0);
                        coo_mat.push(*y-group_sizes[j]+*node_j, *x-group_sizes[i]+node_i, 1.0);
                        degrees[*x-group_sizes[i]+node_i] += 1.0;
                        degrees[*y-group_sizes[j]+*node_j] += 1.0;
                        in_degree[*node_j] -= 1;
                        if i == j {
                            out_degree[node_i] -= 1;
                        }
                        repeats.push((node_i, *node_j));
                    }
                }
            }
        }
        // define ages from partitioning and adjacency matrix as Csr mat
        let mut last_idx = 0;
        let ages: Vec<usize> = partitions  
            .iter()
            .enumerate()
            .flat_map(|(i,x)| {
                let answer = vec![i; *x - last_idx];
                last_idx = *x;
                answer
            })
            .collect();
        let matrix: CsrMatrix<f64> = CsrMatrix::from(&coo_mat);

        NetworkStructure {
            adjacency_matrix: matrix,
            degree: degrees,
            age_brackets: ages
        }
    }

    pub fn new_ba(n: usize, m0: usize, m: usize) -> NetworkStructure {
        //check dimensions correct
        if m0 < m {
            println!("m0 must be greater than m");
        }
        // start with complete sparse graph,
        let mut coo_mat:CooMatrix<f64> = CooMatrix::new(n,n);
        let mut degrees: Vec<f64> = vec![0.0;n];

        // make complete starting graph
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
            degree: degrees,
            age_brackets: Vec::new()
        }
    }

    pub fn new_sbm(n: usize, partitions: Vec<usize>, rates_mat: Vec<Vec<f64>>) -> NetworkStructure {

        // transform rates matrix to probability matrix 
        let prob_mat: Vec<Vec<f64>> = rates_to_probabilities(rates_mat, &partitions); 

        let mut rng: ThreadRng = rand::thread_rng();
        let mut coo_mat: CooMatrix<f64> = CooMatrix::new(n,n);
        let mut degrees: Vec<f64> = vec![0.0;n];

        // loop through partitions and then individuals
        let mut last_idx: [usize;2] = [0,0];
        let mut rand_num: f64;
        
        // loop through lower triangular
        let mut part_i: usize; let mut part_j: usize;
        for i in 0..n {
            for j in 0..i {
                // find which block we are in
                part_i = partitions
                    .iter()
                    .position(|&x| (i/x) < 1)
                    .unwrap();
                part_j = partitions
                    .iter()
                    .position(|&x| (j/x) < 1)
                    .unwrap();
                // randomly generate edges with probability prob_mat
                rand_num = rng.gen();
                if rand_num < prob_mat[part_i][part_j] {
                    coo_mat.push(i, j, 1.0);
                    coo_mat.push(j, i, 1.0);
                    degrees[i] += 1.0;
                    degrees[j] += 1.0;
                }
            }
        }
        
        // define ages from partitioning and adjacency matrix as Csr mat
        last_idx[0] = 0;
        let ages: Vec<usize> = partitions  
            .iter()
            .enumerate()
            .flat_map(|(i,x)| {
                let answer = vec![i; *x - last_idx[0]];
                last_idx[0] = *x;
                answer
            })
            .collect();
        let matrix: CsrMatrix<f64> = CsrMatrix::from(&coo_mat);

        // ensure ages and degree length is consistent
        assert_eq!(ages.len(), degrees.len()); 

        NetworkStructure {
            adjacency_matrix: matrix,
            degree: degrees,
            age_brackets: ages
        }
    }

    pub fn new_sbm_weighted(n: usize, partitions: Vec<usize>, rates_mat: Vec<Vec<f64>>) -> NetworkStructure {
        // unfinished weighting step !!
        // find consecutive group sizes to turn rates to probabilities
        let mut group_sizes: Vec<usize> = partitions
            .windows(2)
            .map(|pair| {
                pair[1] - pair[0]
            })
            .collect();
        group_sizes.insert(0,partitions[0]);
        // transform rates matrix to probability matrix 
        let prob_mat: Vec<Vec<f64>> = rates_mat
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter().map(|rate| {
                    rate / (group_sizes[i] as f64)
                })
                .collect()
            })
            .collect();
        // finish this, do weight calculation here from NB, mu = 4.79, k = 0.54
        let mu: f64 = 4.79; let k: f64 = 0.54;
        let mut rng: ThreadRng = rand::thread_rng();
        let weights: Vec<f64> = NegativeBinomial::new(k, k/(k+mu))
            .unwrap()
            .sample_iter(&mut rng)
            .take(n)
            .map(|x| (x as f64) + 1.0)
            .collect();
        
        let mut coo_mat: CooMatrix<f64> = CooMatrix::new(n,n);
        let mut degrees: Vec<f64> = vec![0.0;n];
        // loop through partitions and then individuals
        let mut last_idx: [usize;2] = [0,0];
        let mut rand_num: f64;
        
        // loop through lower triangular
        let mut part_i: usize; let mut part_j: usize;
        for i in 0..n {
            for j in 0..i {
                // find which block we are in
                part_i = partitions
                    .iter()
                    .position(|&x| (i/x) < 1)
                    .unwrap();
                part_j = partitions
                    .iter()
                    .position(|&x| (j/x) < 1)
                    .unwrap();
                // randomly generate edges with probability prob_mat
                rand_num = rng.gen();
                if rand_num < prob_mat[part_i][part_j] * weights[i] * weights[j] {
                    coo_mat.push(i, j, 1.0);
                    coo_mat.push(j, i, 1.0);
                    degrees[i] += 1.0;
                    degrees[j] += 1.0;
                }
            }
        }
        
        // define ages from partitioning and adjacency matrix as Csr mat
        last_idx[0] = 0;
        let ages: Vec<usize> = partitions  
            .iter()
            .enumerate()
            .flat_map(|(i,x)| {
                let answer = vec![i; *x - last_idx[0]];
                last_idx[0] = *x;
                answer
            })
            .collect();
        let matrix: CsrMatrix<f64> = CsrMatrix::from(&coo_mat);

        // ensure ages and degree length is consistent
        assert_eq!(ages.len(), degrees.len()); 

        NetworkStructure {
            adjacency_matrix: matrix,
            degree: degrees,
            age_brackets: ages
        }
    }
}

impl NetworkProperties {

    pub fn new(network: &NetworkStructure) -> NetworkProperties {
        NetworkProperties { 
            nodal_states: vec![State::Susceptible; network.degree.len()],
            results: Vec::new(),
            result_type: ResultType::SIR,
            outbreak_type: OutbreakType::SIR,
            parameters: vec![0.1,0.2]
        }
    }

    pub fn params(&mut self, params: Vec<f64>) {
        match params.len() {
            2 => {
                self.parameters = params;
                self.outbreak_type = OutbreakType::SIR;
            }
            3 => {
                self.parameters = params;
                self.outbreak_type = OutbreakType::SIRS;
            }
            _ => {
                println!("Enter a vector of 2/3 parameters")
            }
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
        Output { sir: Vec::new(), infections: Vec::new(), network_struct: SerializeableNetwork::new() }
    }
}

#[derive(Debug,Serialize)]
pub struct SerializeableNetwork {
    row_idx: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
    ages: Vec<usize>,
    degrees: Vec<f64>
}

impl SerializeableNetwork {
    
    pub fn new() -> SerializeableNetwork {
        SerializeableNetwork { 
            row_idx: Vec::new(), 
            col_idx: Vec::new(), 
            values: Vec::new(), 
            ages: Vec::new(), 
            degrees: Vec::new() 
        }
    }

    pub fn from(network_structure: &NetworkStructure) -> SerializeableNetwork {
        let coo_mat: CooMatrix<f64> = CooMatrix::from(&network_structure.adjacency_matrix);
        SerializeableNetwork {
            row_idx: coo_mat.row_indices().iter().map(|&x| x).collect(),
            col_idx: coo_mat.col_indices().iter().map(|&x| x).collect(),
            values: coo_mat.values().iter().map(|&x| x).collect(),
            ages: network_structure.age_brackets.clone(),
            degrees: network_structure.degree.clone()
        }
    }
}