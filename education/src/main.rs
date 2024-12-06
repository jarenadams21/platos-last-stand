use rand::{rngs::StdRng, Rng, SeedableRng};
use std::f64::consts::PI;
use ndarray::prelude::*;

/// Constants (representing baseline truths or stable reference points)
const BASELINE: f64 = 1.0;      // A reference "truth level"
const NOISE_LEVEL: f64 = 1e-3;  // Noise magnitude
const LATTICE_SIZE: usize = 20;
const TIME_STEPS: usize = 1000;
const DELTA_T: f64 = 1e-22;      // Time step

/// Probability of stochastic action at a distance: how often a node
/// randomly interacts with a distant node.
const DISTANT_INTERACTION_PROB: f64 = 0.05;

/// Strength of influence when a distant action occurs.
const DISTANT_INFLUENCE_STRENGTH: f64 = 0.1;

/// Probability of a node shifting towards or away from baseline when influenced.
const SHIFT_TOWARDS_BASELINE_PROB: f64 = 0.5;

#[derive(Clone, Copy, Debug)]
struct NodeState {
    /// Represents the node's current "belief intensity"—a scalar field
    /// that might analogize ideological conviction. 
    belief_intensity: f64,
}

impl NodeState {
    fn new(rng: &mut StdRng) -> Self {
        // Initialize with random perturbation around baseline
        let init_val = BASELINE + rng.gen_range(-0.1..0.1);
        NodeState { belief_intensity: init_val }
    }

    /// Apply local random noise as a form of perturbation
    fn apply_noise(&mut self, rng: &mut StdRng) {
        let noise = rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
        self.belief_intensity += noise;
    }

    /// Shift belief intensity towards or away from the baseline
    /// when influenced by a distant node.
    fn distant_influence(&mut self, other_intensity: f64, rng: &mut StdRng) {
        // Check if we push this node towards or away from baseline
        // We can imagine that sometimes distant signals reinforce "truth" and sometimes distort it.
        let direction = if rng.gen::<f64>() < SHIFT_TOWARDS_BASELINE_PROB {
            // Shift towards baseline
            BASELINE
        } else {
            // Shift away from baseline, using the other node's intensity as a reference
            other_intensity
        };

        // Move a fraction of the difference
        let diff = direction - self.belief_intensity;
        self.belief_intensity += diff * DISTANT_INFLUENCE_STRENGTH;
    }

    /// Keep intensities in a reasonable range
    fn clamp(&mut self) {
        // Just to keep numbers manageable
        if self.belief_intensity.is_nan() {
            self.belief_intensity = BASELINE;
        } else {
            self.belief_intensity = self.belief_intensity.max(-10.0).min(10.0);
        }
    }
}

struct Lattice {
    size: usize,
    nodes: Array3<NodeState>,
}

impl Lattice {
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(42);
        let nodes = Array3::from_shape_fn((size, size, size), |_| {
            NodeState::new(&mut rng)
        });

        Lattice { size, nodes }
    }

    fn evolve(&mut self) {
        let mut rng = StdRng::seed_from_u64(1234);

        for _step in 0..TIME_STEPS {
            let nodes_copy = self.nodes.clone();

            // Update nodes in parallel logic:
            // 1. Apply local noise.
            // 2. Occasionally pick a random distant node and influence this node's belief.
            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        let mut node = nodes_copy[[x, y, z]];

                        // Local noise
                        node.apply_noise(&mut rng);

                        // Stochastic action at a distance
                        if rng.gen::<f64>() < DISTANT_INTERACTION_PROB {
                            // Pick a distant node at random
                            let dx = rng.gen_range(0..self.size);
                            let dy = rng.gen_range(0..self.size);
                            let dz = rng.gen_range(0..self.size);

                            let distant_node = nodes_copy[[dx, dy, dz]];
                            node.distant_influence(distant_node.belief_intensity, &mut rng);
                        }

                        // Clamp values
                        node.clamp();

                        self.nodes[[x, y, z]] = node;
                    }
                }
            }

            // One could add global drift or check stabilization after some steps.
        }
    }

    fn measure_divergence(&self) -> f64 {
        // Measure how far on average the lattice is from baseline
        let mut total_diff = 0.0;
        let count = (self.size * self.size * self.size) as f64;
        for node in self.nodes.iter() {
            total_diff += (node.belief_intensity - BASELINE).abs();
        }
        total_diff / count
    }

    fn intensity_variance(&self) -> f64 {
        let mut intensities = Vec::new();
        for node in self.nodes.iter() {
            intensities.push(node.belief_intensity);
        }

        let mean = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let var = intensities.iter().map(|i| (i - mean).powi(2)).sum::<f64>() / (intensities.len() as f64);
        var
    }

    fn print_stats(&self) {
        let divergence = self.measure_divergence();
        let variance = self.intensity_variance();
        println!("Average divergence from baseline: {}", divergence);
        println!("Intensity variance: {}", variance);
    }
}

fn main() {
    let mut lattice = Lattice::new(LATTICE_SIZE);
    println!("Initial state:");
    lattice.print_stats();

    lattice.evolve();
    lattice.evolve();
    lattice.evolve();
    lattice.evolve();
    lattice.evolve();

    println!("Final state after evolution:");
    lattice.print_stats();

    println!("Simulation complete. The system used a stochastic action-at-a-distance approach rather than a Schrödinger representation.");
}