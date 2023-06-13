// use ode_solvers::dopri5::*;
// use ode_solvers::*;
// use crate::random_graphs::*;
// use nalgebra::{DMatrix,DVector};

// type State = DVector<f64>;
// type Time = f64;

// impl ode_solvers::System<State> for Graph {
    
//     fn system(&self, _t: Time, y: &State, dy: &mut State) {
//         let n: usize = y.len()/3;

//         let mut tmp: Vec<DVector<f64>> = vec![
//             - self.adjacency_matrix.clone() * y.column(0).component_mul(&y.column(1))
//         ];
//         tmp.push(
//             self.adjacency_matrix.clone() * y.column(0).component_mul(&y.column(1)) - 
//             self.parameters[0] * y.column(1)
//         );
//         tmp.push(
//             self.adjacency_matrix.clone() * y.column(0).component_mul(&y.column(1)) -
//             self.parameters[0] * y.column(1)
//         );
//         for i in 0..3 {
//             for j in 0..n {
//                 dy[(i*n + j)] = tmp[i][(j)];
//             }
//         }
//     }
// }