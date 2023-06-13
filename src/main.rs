use networks::random_graphs::*;
use networks::odes::*;
use ode_solvers::*;
use ode_solvers::dopri5::*;


type State = DVector<f64>;
type Time = f64;

fn main() {
    // define network
    const n: usize = 10;
    let network = Network::new_ba(n,3,3);
    // println!("{:?}", network.degree);

    // run SIR on network
    network.parameters.push(1.0f64);
    let y0 = State::from_fn(3*n, |r,c| if r < n {0.9} else if r < 2*n {0.1} else {0.0});
    
    let mut stepper = Dopri5::new(network, 0.0, y0, 1.0, 100.0);
}

