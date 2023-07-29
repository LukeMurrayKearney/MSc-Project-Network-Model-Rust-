use rand_distr::{Binomial, Distribution};
use rand::rngs::ThreadRng;

pub fn multinomial_sample(n: usize, ps: Vec<f64>, rng: &mut ThreadRng) -> Vec<usize> {
    let mut x_sum: usize = 0;
    let mut p_sum: f64 = 0.0;
    ps.iter()
        .map(|&p| {
            // make each binomial sample conditional on the last
            let bin = Binomial::new((n - x_sum) as u64, p / (1.0 - p_sum)).unwrap();
            p_sum += p;
            let x: usize = bin.sample(rng) as usize;
            x_sum += x;
            x
        })
        .collect()
}