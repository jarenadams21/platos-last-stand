use ndarray::prelude::*;
use num_complex::Complex;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use rand::Rng;
use std::f64::consts::PI;

/// Physical constants
const HBAR: f64 = 1.0545718e-34; // Reduced Planck constant (J·s)
const MU_B: f64 = 9.274009994e-24; // Bohr magneton (J/T)
const KB: f64 = 1.380649e-23; // Boltzmann constant (J/K)
const GAMMA: f64 = 1.760859e11; // Gyromagnetic ratio (rad·s⁻¹·T⁻¹)
const LATTICE_SIZE: usize = 20; // Lattice dimensions (20x20x20)
const TIME_STEPS: usize = 100;
const DELTA_T: f64 = 1e-22; // Time step (s)
const TEMPERATURE: f64 = 300.0; // Temperature (K)
const EXTERNAL_FIELD: [f64; 3] = [0.0, 0.0, 1.0]; // External magnetic field (T)
const J_EXCHANGE: f64 = 1e-21; // Exchange interaction energy (J)

/// Spinor representing a spin-½ particle
#[derive(Clone, Copy, Debug)]
struct Spinor {
    up: Complex<f64>,
    down: Complex<f64>,
}

impl Spinor {
    /// Creates a new normalized spinor
    fn new(up: Complex<f64>, down: Complex<f64>) -> Self {
        let mut spinor = Spinor { up, down };
        spinor.normalize();
        spinor
    }

    /// Initializes a spinor in a random orientation
    fn random(rng: &mut StdRng) -> Self {
        let theta = rng.gen_range(0.0..PI);
        let phi = rng.gen_range(0.0..2.0 * PI);

        let up = Complex::new((theta / 2.0).cos(), 0.0);
        let down = Complex::new(
            (theta / 2.0).sin() * phi.cos(),
            (theta / 2.0).sin() * phi.sin(),
        );

        Spinor::new(up, down)
    }

    /// Normalizes the spinor
    fn normalize(&mut self) {
        let norm = (self.up.norm_sqr() + self.down.norm_sqr()).sqrt();
        self.up /= norm;
        self.down /= norm;
    }

    /// Returns the spin vector components (expectation values)
    fn spin_vector(&self) -> [f64; 3] {
        let sx = 2.0 * (self.up.conj() * self.down).re;
        let sy = 2.0 * (self.up.conj() * self.down).im;
        let sz = (self.up.conj() * self.up - self.down.conj() * self.down).re;
        [sx, sy, sz]
    }
}

/// Metric tensor representing spacetime curvature
#[derive(Clone, Copy, Debug)]
struct MetricTensor {
    g: [[f64; 3]; 3], // 3x3 metric tensor for simplicity (t, x, y)
}

impl MetricTensor {
    /// Creates a new flat metric tensor (Minkowski spacetime)
    fn new_flat() -> Self {
        let mut g = [[0.0; 3]; 3];
        g[0][0] = 1.0;  // Time component
        g[1][1] = -1.0; // x component
        g[2][2] = -1.0; // y component
        MetricTensor { g }
    }

    /// Introduces perturbations to simulate curvature
    fn perturb(&mut self, rng: &mut StdRng, strength: f64) {
        for i in 0..3 {
            for j in i..3 {
                let perturbation = rng.gen_range(-strength..strength);
                self.g[i][j] += perturbation;
                if i != j {
                    self.g[j][i] = self.g[i][j]; // Ensure symmetry
                }
            }
        }
    }

    /// Computes the inverse metric tensor
    fn inverse(&self) -> [[f64; 3]; 3] {
        let mut inv_g = [[0.0; 3]; 3];
        // Calculate determinant
        let det = self.g[0][0] * (self.g[1][1] * self.g[2][2] - self.g[1][2] * self.g[2][1])
            - self.g[0][1] * (self.g[1][0] * self.g[2][2] - self.g[1][2] * self.g[2][0])
            + self.g[0][2] * (self.g[1][0] * self.g[2][1] - self.g[1][1] * self.g[2][0]);

        let inv_det = 1.0 / det;

        // Compute inverse using adjugate matrix
        inv_g[0][0] = inv_det
            * (self.g[1][1] * self.g[2][2] - self.g[1][2] * self.g[2][1]);
        inv_g[0][1] = inv_det
            * (self.g[0][2] * self.g[2][1] - self.g[0][1] * self.g[2][2]);
        inv_g[0][2] = inv_det
            * (self.g[0][1] * self.g[1][2] - self.g[0][2] * self.g[1][1]);
        inv_g[1][0] = inv_det
            * (self.g[1][2] * self.g[2][0] - self.g[1][0] * self.g[2][2]);
        inv_g[1][1] = inv_det
            * (self.g[0][0] * self.g[2][2] - self.g[0][2] * self.g[2][0]);
        inv_g[1][2] = inv_det
            * (self.g[0][2] * self.g[1][0] - self.g[0][0] * self.g[1][2]);
        inv_g[2][0] = inv_det
            * (self.g[1][0] * self.g[2][1] - self.g[1][1] * self.g[2][0]);
        inv_g[2][1] = inv_det
            * (self.g[0][1] * self.g[2][0] - self.g[0][0] * self.g[2][1]);
        inv_g[2][2] = inv_det
            * (self.g[0][0] * self.g[1][1] - self.g[0][1] * self.g[1][0]);

        inv_g
    }

    /// Computes the Christoffel symbols of the first kind
    fn christoffel_symbols(&self, lattice: &Lattice, x: usize, y: usize, z: usize) -> [[[f64; 3]; 3]; 3] {
        let mut gamma = [[[0.0; 3]; 3]; 3];
        let inv_g = self.inverse();

        for l in 0..3 {
            for m in 0..3 {
                for n in 0..3 {
                    let mut sum = 0.0;
                    for k in 0..3 {
                        let dg_km = self.partial_derivative(k, m, lattice, x, y, z);
                        let dg_kn = self.partial_derivative(k, n, lattice, x, y, z);
                        let dg_mn = self.partial_derivative(m, n, lattice, x, y, z);

                        sum += 0.5 * inv_g[l][k] * (dg_km + dg_kn - dg_mn);
                    }
                    gamma[l][m][n] = sum;
                }
            }
        }
        gamma
    }

    /// Computes the Riemann curvature tensor
    fn riemann_tensor(&self, lattice: &Lattice, x: usize, y: usize, z: usize) -> [[[[f64; 3]; 3]; 3]; 3] {
        let mut riemann = [[[[0.0; 3]; 3]; 3]; 3];
        let gamma = self.christoffel_symbols(lattice, x, y, z);

        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let d_gamma_ilk = self.partial_derivative_gamma(&gamma, i, l, k, j, lattice, x, y, z);
                        let d_gamma_ijk = self.partial_derivative_gamma(&gamma, i, j, k, l, lattice, x, y, z);

                        let mut sum_m = 0.0;
                        for mon in 0..3 {
                            sum_m += gamma[i][j][mon] * gamma[mon][l][k] - gamma[i][l][mon] * gamma[mon][j][k];
                        }

                        riemann[i][j][k][l] = d_gamma_ilk - d_gamma_ijk + sum_m;
                    }
                }
            }
        }
        riemann
    }

    /// Computes the Ricci tensor by contracting the Riemann tensor
    fn ricci_tensor(&self, lattice: &Lattice, x: usize, y: usize, z: usize) -> [[f64; 3]; 3] {
        let mut ricci = [[0.0; 3]; 3];
        let riemann = self.riemann_tensor(lattice, x, y, z);

        for i in 0..3 {
            for k in 0..3 {
                let mut sum_m = 0.0;
                for m in 0..3 {
                    sum_m += riemann[m][i][k][m];
                }
                ricci[i][k] = sum_m;
            }
        }
        ricci
    }

    /// Computes the Ricci scalar by contracting the Ricci tensor
    fn ricci_scalar(&self, lattice: &Lattice, x: usize, y: usize, z: usize) -> f64 {
        let ricci = self.ricci_tensor(lattice, x, y, z);
        let inv_g = self.inverse();
        let mut scalar = 0.0;

        for i in 0..3 {
            for j in 0..3 {
                scalar += inv_g[i][j] * ricci[i][j];
            }
        }
        scalar
    }

    /// Computes the partial derivative of the metric tensor component g_ab with respect to coordinate c
    fn partial_derivative(&self, a: usize, b: usize, lattice: &Lattice, x: usize, y: usize, z: usize) -> f64 {
        let delta = 1; // Grid spacing (assuming unit grid for simplicity)
        let size = lattice.size;

        let prev_index = |i: usize| if i == 0 { size - 1 } else { i - 1 };
        let next_index = |i: usize| if i == size - 1 { 0 } else { i + 1 };

        let mut derivative = 0.0;

        // Finite difference approximation
        match a {
            0 => {
                // Time derivative (we don't vary time in this simulation)
                derivative = 0.0;
            }
            1 => {
                // Spatial derivative with respect to x
                let g_prev = lattice.metric_tensors[[prev_index(x), y, z]].g[b][a];
                let g_next = lattice.metric_tensors[[next_index(x), y, z]].g[b][a];
                derivative = (g_next - g_prev) as f64 / (2.0 * delta as f64);
            }
            2 => {
                // Spatial derivative with respect to y
                let g_prev = lattice.metric_tensors[[x, prev_index(y), z]].g[b][a];
                let g_next = lattice.metric_tensors[[x, next_index(y), z]].g[b][a];
                derivative = (g_next - g_prev) as f64 / (2.0 * delta as f64);
            }
            _ => {}
        }

        derivative
    }

    /// Computes the partial derivative of a Christoffel symbol
    fn partial_derivative_gamma(
        &self,
        gamma: &[[[f64; 3]; 3]; 3],
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        lattice: &Lattice,
        x: usize,
        y: usize,
        z: usize,
    ) -> f64 {
        let delta = 1; // Grid spacing (assuming unit grid for simplicity)
        let size = lattice.size;

        let prev_index = |i: usize| if i == 0 { size - 1 } else { i - 1 };
        let next_index = |i: usize| if i == size - 1 { 0 } else { i + 1 };

        let mut derivative = 0.0;

        // Finite difference approximation
        match l {
            0 => {
                // Time derivative (we don't vary time in this simulation)
                derivative = 0.0;
            }
            1 => {
                // Spatial derivative with respect to x
                let gamma_prev = self.christoffel_symbols(lattice, prev_index(x), y, z)[i][j][k];
                let gamma_next = self.christoffel_symbols(lattice, next_index(x), y, z)[i][j][k];
                derivative = (gamma_next - gamma_prev) / (2.0 * delta as f64);
            }
            2 => {
                // Spatial derivative with respect to y
                let gamma_prev = self.christoffel_symbols(lattice, x, prev_index(y), z)[i][j][k];
                let gamma_next = self.christoffel_symbols(lattice, x, next_index(y), z)[i][j][k];
                derivative = (gamma_next - gamma_prev) / (2.0 * delta as f64);
            }
            _ => {}
        }

        derivative
    }
}

/// Lattice representing the spin system
struct Lattice {
    spins: Array3<Spinor>,
    metric_tensors: Array3<MetricTensor>,
    size: usize,
}

impl Lattice {
    /// Initializes the lattice with random spins and flat metric tensors
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(0);
        let spins = Array3::from_shape_fn((size, size, size), |_| Spinor::random(&mut rng));
        let metric_tensors = Array3::from_shape_fn((size, size, size), |_| {
            let mut metric = MetricTensor::new_flat();
            metric.perturb(&mut rng, 0.01);
            metric
        });
        Lattice {
            spins,
            metric_tensors,
            size,
        }
    }

    /// Evolves the lattice over time
    fn evolve(&mut self) {
        let mut rng = StdRng::seed_from_u64(1);

        for _ in 0..TIME_STEPS {
            let spins_copy = self.spins.clone();

            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        let spin = spins_copy[[x, y, z]];

                        // Compute local effective field
                        let exchange_field =
                            self.compute_exchange_field(x, y, z, &spins_copy);
                        let thermal_field = self.compute_thermal_field(&mut rng);
                        let curvature_effect = self.compute_curvature_effect(x, y, z);

                        let total_field = [
                            EXTERNAL_FIELD[0]
                                + exchange_field[0]
                                + thermal_field[0]
                                + curvature_effect[0],
                            EXTERNAL_FIELD[1]
                                + exchange_field[1]
                                + thermal_field[1]
                                + curvature_effect[1],
                            EXTERNAL_FIELD[2]
                                + exchange_field[2]
                                + thermal_field[2]
                                + curvature_effect[2],
                        ];

                        // Construct Hamiltonian
                        let hamiltonian = self.construct_hamiltonian(&total_field);

                        // Time evolution operator U = exp(-i * H * Δt / ħ)
                        let u_matrix = self.time_evolution_operator(&hamiltonian);

                        // Apply time evolution
                        let new_spinor = self.apply_time_evolution(&spin, &u_matrix);

                        self.spins[[x, y, z]] = new_spinor;
                    }
                }
            }

            // Particle creation and annihilation processes
            self.particle_creation_annihilation(&mut rng);

            // Apply group symmetry operations
            self.apply_symmetry_operations();
        }
    }

    /// Computes the exchange field from neighboring spins
    fn compute_exchange_field(
        &self,
        x: usize,
        y: usize,
        z: usize,
        spins: &Array3<Spinor>,
    ) -> [f64; 3] {
        let mut field = [0.0, 0.0, 0.0];
        let neighbors = self.get_neighbor_indices(x, y, z);

        for &(nx, ny, nz) in &neighbors {
            let neighbor_spin = spins[[nx, ny, nz]];
            let spin_vector = neighbor_spin.spin_vector();
            for i in 0..3 {
                field[i] += spin_vector[i];
            }
        }

        for i in 0..3 {
            field[i] *= J_EXCHANGE / (MU_B * neighbors.len() as f64);
        }

        field
    }

    /// Computes thermal fluctuations
    fn compute_thermal_field(&self, rng: &mut StdRng) -> [f64; 3] {
        let std_dev = (KB * TEMPERATURE / (MU_B)).sqrt();
        let normal = Normal::new(0.0, std_dev).unwrap();
        [
            normal.sample(rng),
            normal.sample(rng),
            normal.sample(rng),
        ]
    }

    /// Computes curvature effects on the local field
    fn compute_curvature_effect(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let metric = self.metric_tensors[[x, y, z]];
        let ricci_scalar = metric.ricci_scalar(&self, x, y, z);

        // Coupling between spin and curvature
        let coupling_constant = 1e-30; // Adjust as appropriate
        let curvature_field = [
            0.0,
            0.0,
            coupling_constant * ricci_scalar,
        ];

        curvature_field
    }

    /// Constructs the Hamiltonian matrix for a spin
    fn construct_hamiltonian(&self, field: &[f64; 3]) -> [[Complex<f64>; 2]; 2] {
        let gamma = GAMMA;
        let hx = field[0];
        let hy = field[1];
        let hz = field[2];

        let h11 = -gamma * hz / 2.0;
        let h22 = gamma * hz / 2.0;

        let hx_complex = Complex::new(hx, 0.0);
        let hy_complex = Complex::new(hy, 0.0);

        let h12 = -gamma * (hx_complex - Complex::<f64>::i() * hy_complex) / 2.0;
        let h21 = -gamma * (hx_complex + Complex::<f64>::i() * hy_complex) / 2.0;

        let hamiltonian = [
            [Complex::new(h11, 0.0), h12],
            [h21, Complex::new(h22, 0.0)],
        ];

        hamiltonian
    }

    /// Computes the time evolution operator U = exp(-i * H * Δt / ħ)
    fn time_evolution_operator(
        &self,
        hamiltonian: &[[Complex<f64>; 2]; 2],
    ) -> [[Complex<f64>; 2]; 2] {
        let factor = Complex::new(0.0, -DELTA_T / HBAR);
        let h_scaled = [
            [
                hamiltonian[0][0] * factor,
                hamiltonian[0][1] * factor,
            ],
            [
                hamiltonian[1][0] * factor,
                hamiltonian[1][1] * factor,
            ],
        ];
        matrix_exponential(&h_scaled)
    }

    /// Applies the time evolution to a spinor
    fn apply_time_evolution(
        &self,
        spinor: &Spinor,
        u_matrix: &[[Complex<f64>; 2]; 2],
    ) -> Spinor {
        let new_up = u_matrix[0][0] * spinor.up + u_matrix[0][1] * spinor.down;
        let new_down = u_matrix[1][0] * spinor.up + u_matrix[1][1] * spinor.down;
        Spinor::new(new_up, new_down)
    }

    /// Particle creation and annihilation processes
    fn particle_creation_annihilation(&mut self, rng: &mut StdRng) {
        // For simplicity, we'll model a basic stochastic creation and annihilation process
        let creation_probability = 0.01;
        let annihilation_probability = 0.01;

        for spin in self.spins.iter_mut() {
            if rng.gen_bool(creation_probability) {
                // Particle creation (random spin state)
                *spin = Spinor::random(rng);
            } else if rng.gen_bool(annihilation_probability) {
                // Particle annihilation (spin state becomes zero)
                *spin = Spinor::new(Complex::new(0.0, 0.0), Complex::new(0.0, 0.0));
            }
        }
    }

    /// Applies symmetry operations based on finite groups
    fn apply_symmetry_operations(&mut self) {
        // Placeholder for group theory operations
        // In a full implementation, apply group elements to the lattice
    }

    /// Retrieves the indices of neighboring spins (with periodic boundary conditions)
    fn get_neighbor_indices(&self, x: usize, y: usize, z: usize) -> Vec<(usize, usize, usize)> {
        let size = self.size;
        vec![
            ((x + size - 1) % size, y, z),
            ((x + 1) % size, y, z),
            (x, (y + size - 1) % size, z),
            (x, (y + 1) % size, z),
            (x, y, (z + size - 1) % size),
            (x, y, (z + 1) % size),
        ]
    }

    /// Calculates the average magnetization of the lattice
    fn calculate_magnetization(&self) -> [f64; 3] {
        let mut total = [0.0, 0.0, 0.0];
        let num_spins = self.size.pow(3) as f64;

        for spin in self.spins.iter() {
            let spin_vector = spin.spin_vector();
            for i in 0..3 {
                total[i] += spin_vector[i];
            }
        }

        for i in 0..3 {
            total[i] /= num_spins;
        }

        total
    }
}

/// Exponential of a 2x2 matrix using scaling and squaring with Padé approximation
fn matrix_exponential(a: &[[Complex<f64>; 2]; 2]) -> [[Complex<f64>; 2]; 2] {
    // Implemented based on expm function in numerical libraries
    // For 2x2 matrices, we can compute the exponential exactly

    let a00 = a[0][0];
    let a01 = a[0][1];
    let a10 = a[1][0];
    let a11 = a[1][1];

    let trace = a00 + a11;
    let delta = (a00 - a11).powi(2) + 4.0 * a01 * a10;
    let sqrt_delta = delta.sqrt();

    let exp_half_trace = (trace / 2.0).exp();

    let cosh = (sqrt_delta / 2.0).cosh();
    let sinh = (sqrt_delta / 2.0).sinh();

    let factor = if sqrt_delta != Complex::new(0.0, 0.0) {
        sinh / sqrt_delta
    } else {
        Complex::new(0.5, 0.0)
    };

    let exp_a = [
        [
            exp_half_trace * (cosh + factor * (a00 - a11) / 2.0),
            exp_half_trace * factor * a01 * 2.0,
        ],
        [
            exp_half_trace * factor * a10 * 2.0,
            exp_half_trace * (cosh - factor * (a00 - a11) / 2.0),
        ],
    ];

    exp_a
}

fn main() {
    // Initialize the lattice
    let mut lattice = Lattice::new(LATTICE_SIZE);

    // Initial magnetization
    let initial_magnetization = lattice.calculate_magnetization();
    println!("Initial Magnetization: {:?}", initial_magnetization);

    // Evolve the lattice
    lattice.evolve();

    // Final magnetization
    let final_magnetization = lattice.calculate_magnetization();
    println!("Final Magnetization: {:?}", final_magnetization);
}