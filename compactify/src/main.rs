use piston_window::*;
use nalgebra::{Complex, DMatrix, DVector};
use rand::Rng;
use std::collections::HashSet;
use ndarray::{Array1, Array2};
use ndarray_linalg::{Eig, Eigh, UPLO};
use num_traits::Zero;
use nalgebra::ComplexField;

/// Represents a single particle with position, velocity, and internal quantum state.
#[derive(Clone)]
struct Particle {
    position: [f64; 2],
    velocity: [f64; 2],
    internal_state: DVector<Complex<f64>>,
    entangled_partner: Option<usize>, // Index of the entangled partner particle
    charged: bool,                    // Indicates if the particle is charged
}

impl Particle {
    /// Creates a new particle at a given position with a random initial internal state.
    fn new(x: f64, y: f64, state_dimension: usize) -> Self {
        // Initialize internal state as a normalized random complex vector
        let mut rng = rand::thread_rng();
        let mut state = DVector::from_fn(state_dimension, |_, _| {
            Complex::new(rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5)
        });

        // Normalize the state vector
        let norm = state.norm();
        state /= Complex::from(norm);

        Particle {
            position: [x, y],
            velocity: [0.0, 0.0],
            internal_state: state,
            entangled_partner: None,
            charged: false, // Initialize as not charged
        }
    }

    /// Updates the particle's position based on its velocity.
    fn update_position(&mut self, dt: f64, width: f64, height: f64) {
        self.position[0] += self.velocity[0] * dt;
        self.position[1] += self.velocity[1] * dt;

        // Keep particles within bounds (wrap around)
        if self.position[0] < 0.0 {
            self.position[0] += width;
        } else if self.position[0] > width {
            self.position[0] -= width;
        }

        if self.position[1] < 0.0 {
            self.position[1] += height;
        } else if self.position[1] > height {
            self.position[1] -= height;
        }
    }

    /// Gets a color representation of the internal state for visualization.
    fn get_color(&self) -> [f32; 4] {
        if self.charged {
            // Charged particles are blue
            [0.0, 0.0, 1.0, 1.0]
        } else {
            // Map the probability amplitudes to color components
            let prob = self.internal_state.map(|c| c.norm_sqr());

            // Assume first three components map to RGB
            let r = prob.get(0).cloned().unwrap_or(0.0) as f32;
            let g = prob.get(1).cloned().unwrap_or(0.0) as f32;
            let b = prob.get(2).cloned().unwrap_or(0.0) as f32;

            // Normalize to [0,1]
            let total = r + g + b;
            let r = if total > 0.0 { r / total } else { 0.0 };
            let g = if total > 0.0 { g / total } else { 0.0 };
            let b = if total > 0.0 { b / total } else { 0.0 };

            [r, g, b, 1.0] // Alpha is 1.0 (opaque)
        }
    }
}

/// Represents the simulation environment.
struct Simulation {
    particles: Vec<Particle>,
    width: f64,
    height: f64,
    hamiltonian_individual: DMatrix<Complex<f64>>,
    hamiltonian_joint: DMatrix<Complex<f64>>,
    state_dimension: usize,
}

impl Simulation {
    /// Creates a new simulation with particles initialized in random positions.
    fn new(num_particles: usize, width: f64, height: f64) -> Self {
        let state_dimension = 4; // For SU(4), the dimension is 4
        let mut particles = Vec::with_capacity(num_particles);
        let mut rng = rand::thread_rng();

        for _ in 0..num_particles {
            let x = rng.gen_range(0.0..width);
            let y = rng.gen_range(0.0..height);
            particles.push(Particle::new(x, y, state_dimension));
        }

        // Entangle random pairs of particles
        for i in (0..num_particles).step_by(2) {
            if i + 1 < num_particles {
                particles[i].entangled_partner = Some(i + 1);
                particles[i + 1].entangled_partner = Some(i);
                // Initialize entangled Bell state
                let norm = 1.0 / (2.0f64).sqrt();
                let bell_state = DVector::from_vec(vec![
                    Complex::new(norm, 0.0),
                    Complex::new(0.0, 0.0),
                    Complex::new(0.0, 0.0),
                    Complex::new(norm, 0.0),
                ]);

                particles[i].internal_state = bell_state.clone();
                particles[i + 1].internal_state = bell_state;
            }
        }

        // Define the individual and joint Hamiltonians
        let hamiltonian_individual = Self::generate_su4_hamiltonian();
        let hamiltonian_joint = Self::generate_joint_hamiltonian(&hamiltonian_individual);

        Simulation {
            particles,
            width,
            height,
            hamiltonian_individual,
            hamiltonian_joint,
            state_dimension,
        }
    }

    /// Generates a complex SU(4) Hamiltonian incorporating rotational and magnetic effects.
    fn generate_su4_hamiltonian() -> DMatrix<Complex<f64>> {
        let dim = 4;

        // Define angular frequency for rotation
        let omega = 1.0;

        // Define magnetic field strength
        let b_field = 1.0;

        // Rotational term
        let rotational_term = DMatrix::from_fn(dim, dim, |i, j| {
            if i == j {
                Complex::new(omega * (i as f64), 0.0)
            } else {
                Complex::new(0.0, 0.0)
            }
        });

        // Magnetic term
        let magnetic_term = DMatrix::from_fn(dim, dim, |i, j| {
            if i == j {
                Complex::new(b_field * (i as f64), 0.0)
            } else {
                Complex::new(0.0, 0.0)
            }
        });

        // Total Hamiltonian
        let mut h = rotational_term + magnetic_term;

        // Make it Hermitian
        let h_dagger = h.adjoint();
        h = (&h + &h_dagger) * Complex::from(0.5);

        h
    }

    /// Generates the joint Hamiltonian for entangled particles.
    fn generate_joint_hamiltonian(
        hamiltonian_individual: &DMatrix<Complex<f64>>,
    ) -> DMatrix<Complex<f64>> {
        let identity = DMatrix::<Complex<f64>>::identity(
            hamiltonian_individual.nrows(),
            hamiltonian_individual.ncols(),
        );
        let hamiltonian_joint = kronecker_product(hamiltonian_individual, &identity)
            + kronecker_product(&identity, hamiltonian_individual);
        hamiltonian_joint
    }

    /// Applies an electroweak magnetic pulse to target particles.
    fn apply_pulse(&mut self) {
        let width = self.width;
        for particle in &mut self.particles {
            if particle.position[0] < width / 2.0 {
                // Modify the particle's Hamiltonian or internal state to simulate charging
                particle.charged = true;
            }
        }
    }

    /// Updates the simulation state.
    fn update(&mut self, dt: f64) {
        // Update positions
        for particle in &mut self.particles {
            particle.update_position(dt, self.width, self.height);
        }

        // Handle entangled pairs
        let mut processed = HashSet::new();
        let len = self.particles.len();
        for i in 0..len {
            if let Some(partner_index) = self.particles[i].entangled_partner {
                if !processed.contains(&i) {
                    // Evolve joint state
                    let particle_a = &self.particles[i];
                    let particle_b = &self.particles[partner_index];

                    let joint_state = kronecker_state(
                        &particle_a.internal_state,
                        &particle_b.internal_state,
                    );

                    // Evolve the joint state
                    let evolved_joint_state =
                        evolve_state(&joint_state, &self.hamiltonian_joint, dt);

                    // Update particles with the new joint state
                    let (new_state_a, new_state_b) =
                        split_joint_state(&evolved_joint_state, self.state_dimension);

                    // Update particles
                    self.particles[i].internal_state = new_state_a;
                    self.particles[partner_index].internal_state = new_state_b;

                    processed.insert(i);
                    processed.insert(partner_index);
                }
            }
        }

        // Update unentangled particles
        for i in 0..len {
            if self.particles[i].entangled_partner.is_none() || !processed.contains(&i) {
                let particle = &mut self.particles[i];
                particle.internal_state =
                    evolve_state(&particle.internal_state, &self.hamiltonian_individual, dt);
            }
        }
    }

    /// Renders the simulation onto the window.
    fn render<G: Graphics>(&self, c: &Context, g: &mut G) {
        for particle in &self.particles {
            let color = particle.get_color();
            ellipse(
                color,
                [
                    particle.position[0] - 3.0,
                    particle.position[1] - 3.0,
                    6.0,
                    6.0,
                ],
                c.transform,
                g,
            );
        }
    }
}

/// Computes the Kronecker product of two matrices.
fn kronecker_product(
    a: &DMatrix<Complex<f64>>,
    b: &DMatrix<Complex<f64>>,
) -> DMatrix<Complex<f64>> {
    let rows = a.nrows() * b.nrows();
    let cols = a.ncols() * b.ncols();
    let mut result = DMatrix::zeros(rows, cols);

    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            let a_elem = a[(i, j)];
            for k in 0..b.nrows() {
                for l in 0..b.ncols() {
                    result[(i * b.nrows() + k, j * b.ncols() + l)] = a_elem * b[(k, l)];
                }
            }
        }
    }
    result
}

/// Computes the Kronecker product of two state vectors.
fn kronecker_state(
    a: &DVector<Complex<f64>>,
    b: &DVector<Complex<f64>>,
) -> DVector<Complex<f64>> {
    let len = a.len() * b.len();
    let mut result = DVector::zeros(len);

    for i in 0..a.len() {
        for j in 0..b.len() {
            result[i * b.len() + j] = a[i] * b[j];
        }
    }
    result
}

fn split_joint_state(
    joint_state: &DVector<Complex<f64>>,
    state_dimension: usize,
) -> (DVector<Complex<f64>>, DVector<Complex<f64>>) {
    // Convert joint_state to ndarray
    let joint_state_array = Array1::from_shape_vec(
        joint_state.len(),
        joint_state.iter().cloned().collect(),
    )
    .expect("Failed to create ndarray from joint_state");

    // Construct the joint density matrix rho_joint = |psi><psi|
    let psi = joint_state_array
        .clone()
        .into_shape((joint_state.len(), 1))
        .unwrap();
    let rho_joint = psi.dot(&psi.mapv(|x| x.conj()).t());

    // Partial trace over the second subsystem
    let dim = state_dimension;
    let mut rho_a = Array2::<Complex<f64>>::zeros((dim, dim));
    let mut rho_b = Array2::<Complex<f64>>::zeros((dim, dim));

    for i in 0..dim {
        for j in 0..dim {
            let mut sum_a = Complex::zero();
            let mut sum_b = Complex::zero();
            for k in 0..dim {
                sum_a += rho_joint[[i * dim + k, j * dim + k]];
                sum_b += rho_joint[[k * dim + i, k * dim + j]];
            }
            rho_a[[i, j]] = sum_a;
            rho_b[[i, j]] = sum_b;
        }
    }

    // Ensure rho_a and rho_b are Hermitian
    rho_a = (&rho_a + &rho_a.t().mapv(|x| x.conj())) * Complex::from(0.5);
    rho_b = (&rho_b + &rho_b.t().mapv(|x| x.conj())) * Complex::from(0.5);

    // Compute eigenvalues and eigenvectors using Eigh
    let (eigenvalues_a, eigenvectors_a) = rho_a
        .eigh(UPLO::Lower)
        .expect("Eigenvalue decomposition failed");
    let (eigenvalues_b, eigenvectors_b) = rho_b
        .eigh(UPLO::Lower)
        .expect("Eigenvalue decomposition failed");

    // Find the eigenvector corresponding to the largest eigenvalue
    let max_idx_a = eigenvalues_a
        .iter()
        .enumerate()
        .max_by(|(_, val_a), (_, val_b)| {
            val_a.norm1().partial_cmp(&val_b.norm1()).unwrap()
        })
        .unwrap()
        .0;
    let new_state_a_array = eigenvectors_a.column(max_idx_a).to_owned();

    let max_idx_b = eigenvalues_b
        .iter()
        .enumerate()
        .max_by(|(_, val_a), (_, val_b)| {
            val_a.norm1().partial_cmp(&val_b.norm1()).unwrap()
        })
        .unwrap()
        .0;
    let new_state_b_array = eigenvectors_b.column(max_idx_b).to_owned();

    // Convert back to nalgebra vectors
    let new_state_a = DVector::from_vec(new_state_a_array.to_vec());
    let new_state_b = DVector::from_vec(new_state_b_array.to_vec());

    // Normalize the state vectors
    let norm_a = new_state_a.norm();
    let norm_b = new_state_b.norm();

    (
        new_state_a / Complex::from(norm_a),
        new_state_b / Complex::from(norm_b),
    )
}

fn evolve_state(
    state: &DVector<Complex<f64>>,
    hamiltonian: &DMatrix<Complex<f64>>,
    dt: f64,
) -> DVector<Complex<f64>> {
    // Convert nalgebra structures to ndarray
    let hamiltonian_array = Array2::from_shape_vec(
        (hamiltonian.nrows(), hamiltonian.ncols()),
        hamiltonian.iter().cloned().collect(),
    )
    .expect("Failed to create ndarray from hamiltonian");

    // Compute eigenvalues and eigenvectors
    let (eigenvalues, eigenvectors) = hamiltonian_array
        .eigh(UPLO::Lower)
        .expect("Eigenvalue decomposition failed");

    // Compute the diagonal matrix of exponentials
    let exp_eigenvalues = Array1::from_shape_fn(eigenvalues.len(), |i| {
        (-Complex::i() * eigenvalues[i] * dt).exp()
    });

    // Construct the exponential matrix U = V * exp(D) * V^H
    let exp_diag = Array2::from_diag(&exp_eigenvalues);
    let u = eigenvectors.dot(&exp_diag).dot(&eigenvectors.t().mapv(|x| x.conj()));

    // Convert the state vector
    let state_array = Array1::from_shape_vec(
        state.len(),
        state.iter().cloned().collect(),
    )
    .expect("Failed to create ndarray from state vector");

    // Evolve the state
    let new_state_array = u.dot(&state_array);

    // Convert back to nalgebra vector
    let new_state = DVector::from_vec(new_state_array.to_vec());

    // Normalize the state vector
    let norm = new_state.norm();
    new_state / Complex::from(norm)
}

fn main() {
    // Window dimensions
    let (width, height) = (800.0, 600.0);

    // Create a new Piston window
    let mut window: PistonWindow =
        WindowSettings::new(
            "Quantum SU(4) Simulation with Entanglement",
            [width as u32, height as u32],
        )
        .exit_on_esc(true)
        .build()
        .expect("Failed to create window.");

    // Simulation parameters
    let num_particles = 5; // Reduced for performance
    let dt = 1.0;

    // Create the simulation
    let mut simulation = Simulation::new(num_particles, width, height);

    // Apply the pulse at the beginning
    simulation.apply_pulse();

    // Event loop
    while let Some(event) = window.next() {
        // Update the simulation
        simulation.update(dt);
        // simulation.apply_pulse();

        // Render the simulation
        window.draw_2d(&event, |c, g, _| {
            // Clear the screen to black
            clear([0.0; 4], g);
            simulation.render(&c, g);
        });
    }
}