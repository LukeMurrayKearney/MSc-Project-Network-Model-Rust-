extern crate nalgebra as na;
use na::{DMatrix, DVector};
extern crate random_choice;
use self::random_choice::random_choice;

#[derive(Debug)]
pub struct Network {
    pub adjacency_matrix: DMatrix<f64>,
    pub degree: DVector<f64>,
    pub parameters: DVector<f64>,
}

impl Network {
    pub fn new_ba(n: usize, m0: usize, m: usize) -> Network {
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
            parameters: DVector::zeros(0)
        }
    }
}
