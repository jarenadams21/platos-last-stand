// Imports
use num_complex::Complex;
use ndarray::prelude::*;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::Normal;
use plotters::prelude::*;
use kiss3d::window::Window;
use kiss3d::light::Light;
use nalgebra::{Point3, Vector3};

// Constants
const HBAR: f64 = 1.0545718e-34;        // Reduced Planck constant in J·s
const MU_B: f64 = 9.274009994e-24;      // Bohr magneton in J/T
const KB: f64 = 1.380649e-23;           // Boltzmann constant in J/K
const J_EXCHANGE: f64 = 1e-21;          // Exchange interaction energy in J
const LATTICE_SIZE: usize = 20;         // Size of the lattice (20x20x20)
const TIME_STEPS: usize = 1000;         // Number of time steps
const DELTA_T: f64 = 1e-22;             // Time step in seconds
const TEMPERATURE: f64 = 1e2;           // Temperature in Kelvin
const EXTERNAL_FIELD: f64 = 0.0;        // External magnetic field in Tesla
const CMB_TEMPERATURE: f64 = 2.725;     // Cosmic Microwave Background temperature in Kelvin
const LATTICE_SPACING: f64 = 1e-10;     // Lattice spacing in meters (e.g., 0.1 nm)

// Spin struct representing a quantum spin state using density matrices
#[derive(Clone, Copy, Debug)]
struct Spin {
    density_matrix: [[Complex<f64>; 2]; 2],
}

impl Spin {
    // Pauli matrices
    const PAULI_X: [[Complex<f64>; 2]; 2] = [
        [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
        [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
    ];
    const PAULI_Y: [[Complex<f64>; 2]; 2] = [
        [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
        [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)],
    ];
    const PAULI_Z: [[Complex<f64>; 2]; 2] = [
        [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)],
    ];

    /// Initialize a spin with a random pure state density matrix
    fn new_random(rng: &mut StdRng) -> Self {
        // Random angles
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let phi = rng.gen_range(0.0..(2.0 * std::f64::consts::PI));

        // Bloch sphere representation
        let nx = theta.sin() * phi.cos();
        let ny = theta.sin() * phi.sin();
        let nz = theta.cos();

        // Convert to density matrix
        let density_matrix = [
            [
                Complex::new((1.0 + nz) / 2.0, 0.0),
                Complex::new((nx - Complex::<f64>::i() * ny) / 2.0, 0.0),
            ],
            [
                Complex::new((nx + Complex::<f64>::i() * ny) / 2.0, 0.0),
                Complex::new((1.0 - nz) / 2.0, 0.0),
            ],
        ];

        Spin { density_matrix }
    }

    /// Calculate expectation value for an operator
    fn expectation(&self, operator: &[[Complex<f64>; 2]; 2]) -> f64 {
        let mut result = Complex::new(0.0, 0.0);
        for i in 0..2 {
            for j in 0..2 {
                result += self.density_matrix[i][j] * operator[j][i];
            }
        }
        result.re
    }

    fn expectation_sx(&self) -> f64 {
        self.expectation(&Self::PAULI_X)
    }

    fn expectation_sy(&self) -> f64 {
        self.expectation(&Self::PAULI_Y)
    }

    fn expectation_sz(&self) -> f64 {
        self.expectation(&Self::PAULI_Z)
    }
}

// Photon struct representing an observable photon
struct Photon {
    position: [f64; 3],
    momentum: [f64; 3],
    polarization: [f64; 3],
    frequency: f64,
}

impl Photon {
    fn new(position: [f64; 3], momentum: [f64; 3], polarization: [f64; 3]) -> Self {
        let frequency = (momentum[0].powi(2) + momentum[1].powi(2) + momentum[2].powi(2)).sqrt() / (2.0 * std::f64::consts::PI);
        Photon {
            position,
            momentum,
            polarization,
            frequency,
        }
    }

    /// Move the photon forward in time
    fn propagate(&mut self, delta_t: f64) {
        for i in 0..3 {
            self.position[i] += self.momentum[i] * delta_t;
        }
    }

    /// Check if the photon is within the lattice bounds
    fn is_in_lattice(&self, size: usize) -> bool {
        self.position.iter().all(|&x| x >= 0.0 && x < size as f64)
    }
}

// Lattice struct representing the 3D lattice of spins
struct Lattice {
    spins: Array3<Spin>,
    size: usize,
    photon: Option<Photon>,
}

impl Lattice {
    /// Initialize a new lattice with spins in random orientations
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(0);
        let spins = Array3::from_shape_fn((size, size, size), |_| {
            Spin::new_random(&mut rng)
        });
        Lattice {
            spins,
            size,
            photon: Some(Photon::new(
                [0.0, size as f64 / 2.0, size as f64 / 2.0], // Starting position at one edge
                [1.0, 0.0, 0.0],   // Momentum along x-axis
                [0.0, 1.0, 0.0],   // Polarization along y-axis
            )),
        }
    }

    /// Simulate one evolution step of the lattice
    fn evolve_step(&mut self) {
        let mut rng = StdRng::seed_from_u64(0);

        let spins_copy = self.spins.clone();
        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    let spin = spins_copy[[x, y, z]];
                    let neighbors = self.get_neighbors(x, y, z);

                    // Compute exchange field from neighbors
                    let mut exchange_field = [0.0, 0.0, 0.0];
                    for neighbor_spin in &neighbors {
                        exchange_field[0] += neighbor_spin.expectation_sx();
                        exchange_field[1] += neighbor_spin.expectation_sy();
                        exchange_field[2] += neighbor_spin.expectation_sz();
                    }

                    // Normalize exchange field and convert to Tesla
                    let num_neighbors = neighbors.len() as f64;
                    exchange_field[0] *= J_EXCHANGE / (MU_B * num_neighbors);
                    exchange_field[1] *= J_EXCHANGE / (MU_B * num_neighbors);
                    exchange_field[2] *= J_EXCHANGE / (MU_B * num_neighbors);

                    // Thermal fluctuations due to lattice temperature (in Tesla)
                    let thermal_std = (2.0 * KB * TEMPERATURE / (MU_B)).sqrt();
                    let normal_dist = Normal::new(0.0, thermal_std).unwrap();
                    let thermal_field = [
                        rng.sample(normal_dist),
                        rng.sample(normal_dist),
                        rng.sample(normal_dist),
                    ];

                    // CMB field fluctuations (in Tesla)
                    let cmb_field = self.calculate_cmb_field(&mut rng);

                    // External magnetic field
                    let external_field = [0.0, 0.0, EXTERNAL_FIELD];

                    // Total effective magnetic field (in Tesla)
                    let total_field = [
                        exchange_field[0] + thermal_field[0] + cmb_field[0] + external_field[0],
                        exchange_field[1] + thermal_field[1] + cmb_field[1] + external_field[1],
                        exchange_field[2] + thermal_field[2] + cmb_field[2] + external_field[2],
                    ];

                    // Magnetic moment of an electron spin (Bohr magneton)
                    const MU_S: f64 = MU_B;

                    // Hamiltonian matrix elements (in Joules)
                    let h11 = -0.5 * MU_S * total_field[2];
                    let h12 = -0.5 * MU_S * (total_field[0] - Complex::<f64>::i() * total_field[1]);
                    let h21 = -0.5 * MU_S * (total_field[0] + Complex::<f64>::i() * total_field[1]);
                    let h22 = 0.5 * MU_S * total_field[2];

                    // Time evolution operator: U = exp(-i * H * Δt / ħ)
                    let exponent = [
                        [Complex::new(h11, 0.0), h12],
                        [h21, Complex::new(h22, 0.0)],
                    ];
                    let exponent =
                        matrix_scalar_multiply(&exponent, Complex::new(0.0, -DELTA_T / HBAR));

                    // Exponentiate the Hamiltonian matrix using Padé approximant
                    let u_matrix = matrix_exponential(&exponent);

                    // Update the density matrix: ρ' = U ρ U†
                    let mut new_density = [[Complex::new(0.0, 0.0); 2]; 2];
                    for i in 0..2 {
                        for j in 0..2 {
                            for k in 0..2 {
                                for l in 0..2 {
                                    new_density[i][j] +=
                                        u_matrix[i][k] * spin.density_matrix[k][l] * u_matrix[j][l].conj();
                                }
                            }
                        }
                    }

                    self.spins[[x, y, z]] = Spin {
                        density_matrix: new_density,
                    };
                }
            }
        }

        // Photon interaction
        if let Some(photon) = &mut self.photon {
            photon.propagate(DELTA_T);

            // Check for interaction with spins
            if photon.is_in_lattice(self.size) {
                let x = photon.position[0].floor() as usize;
                let y = photon.position[1].floor() as usize;
                let z = photon.position[2].floor() as usize;

                // Apply interaction
                let spin = &mut self.spins[[x, y, z]];

                // Model photon-induced spin flip with some probability
                let flip_probability = 0.1; // Adjust as needed
                if rng.gen_bool(flip_probability) {
                    // Flip the spin state
                    let temp = spin.density_matrix[0][0];
                    spin.density_matrix[0][0] = spin.density_matrix[1][1];
                    spin.density_matrix[1][1] = temp;
                    spin.density_matrix[0][1] = -spin.density_matrix[0][1];
                    spin.density_matrix[1][0] = -spin.density_matrix[1][0];
                }
            } else {
                // Remove the photon if it exits the lattice
                self.photon = None;
            }
        }
    }

    /// Calculate the CMB field fluctuations for a spin
    fn calculate_cmb_field(&self, rng: &mut StdRng) -> [f64; 3] {
        // The CMB photons interact weakly but are modeled as an additional thermal field
        let cmb_std = (2.0 * KB * CMB_TEMPERATURE / (MU_B)).sqrt();
        let normal_dist = Normal::new(0.0, cmb_std).unwrap();
        [
            rng.sample(normal_dist),
            rng.sample(normal_dist),
            rng.sample(normal_dist),
        ]
    }

    /// Get the neighboring spins for a given position
    fn get_neighbors(&self, x: usize, y: usize, z: usize) -> Vec<Spin> {
        let mut neighbors = Vec::new();
        let size = self.size;
        let positions = [
            ((x + size - 1) % size, y, z),
            ((x + 1) % size, y, z),
            (x, (y + size - 1) % size, z),
            (x, (y + 1) % size, z),
            (x, y, (z + size - 1) % size),
            (x, y, (z + 1) % size),
        ];
        for &(nx, ny, nz) in &positions {
            neighbors.push(self.spins[[nx, ny, nz]]);
        }
        neighbors
    }
}

/// Matrix multiplication for 2x2 complex matrices
fn matrix_multiply(
    a: &[[Complex<f64>; 2]; 2],
    b: &[[Complex<f64>; 2]; 2],
) -> [[Complex<f64>; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// Matrix exponential using Padé approximant for 2x2 matrices
fn matrix_exponential(a: &[[Complex<f64>; 2]; 2]) -> [[Complex<f64>; 2]; 2] {
    // For small matrices, Padé approximant provides good numerical stability
    let n = 6; // Order of the approximant

    let mut a_power = [
        [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
        [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
    ]; // a^0 = I

    let mut numerator = [
        [Complex::new(0.0, 0.0); 2],
        [Complex::new(0.0, 0.0); 2],
    ];
    let mut denominator = numerator;

    for k in 0..=n {
        if k > 0 {
            a_power = matrix_multiply(&a_power, a);
        }
        let coeff = factorial(2 * n - k) * factorial(k);
        let coeff = Complex::new(coeff as f64, 0.0);

        let term = matrix_scalar_multiply(&a_power, coeff);

        if k % 2 == 0 {
            numerator = matrix_add(&numerator, &term);
            denominator = matrix_add(&denominator, &term);
        } else {
            numerator = matrix_add(&numerator, &matrix_scalar_multiply(&term, Complex::new(-1.0, 0.0)));
            denominator = matrix_add(&denominator, &matrix_scalar_multiply(&term, Complex::new(-1.0, 0.0)));
        }
    }

    // Inverse of denominator
    let det = denominator[0][0] * denominator[1][1] - denominator[0][1] * denominator[1][0];
    let inv_det = det.inv();
    let inverse_denominator = [
        [denominator[1][1] * inv_det, -denominator[0][1] * inv_det],
        [-denominator[1][0] * inv_det, denominator[0][0] * inv_det],
    ];

    matrix_multiply(&inverse_denominator, &numerator)
}

/// Matrix addition for 2x2 complex matrices
fn matrix_add(
    a: &[[Complex<f64>; 2]; 2],
    b: &[[Complex<f64>; 2]; 2],
) -> [[Complex<f64>; 2]; 2] {
    [
        [a[0][0] + b[0][0], a[0][1] + b[0][1]],
        [a[1][0] + b[1][0], a[1][1] + b[1][1]],
    ]
}

/// Scalar multiplication of a 2x2 complex matrix
fn matrix_scalar_multiply(
    a: &[[Complex<f64>; 2]; 2],
    scalar: Complex<f64>,
) -> [[Complex<f64>; 2]; 2] {
    [
        [a[0][0] * scalar, a[0][1] * scalar],
        [a[1][0] * scalar, a[1][1] * scalar],
    ]
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

fn main() {
    // Initialize the lattice
    let mut lattice = Lattice::new(LATTICE_SIZE);

    // Set up the 3D visualization window
    let mut window = Window::new("3D Spin Lattice Visualization");
    window.set_light(Light::StickToCamera);

    // Create a vector to hold the cubes representing spins
    let mut cubes = Vec::new();

    // Add cubes to the window representing each spin
    for x in 0..LATTICE_SIZE {
        for y in 0..LATTICE_SIZE {
            for z in 0..LATTICE_SIZE {
                let mut cube = window.add_cube(0.8, 0.8, 0.8);
                cube.set_local_translation(
                    Vector3::new(x as f32, y as f32, z as f32).into(),
                );
                cubes.push(cube);
            }
        }
    }

    // Main loop
    while window.render() {
        // Evolve the lattice
        lattice.evolve_step();

        // Update the cubes' colors based on spin orientation
        for (idx, cube) in cubes.iter_mut().enumerate() {
            let x = idx / (LATTICE_SIZE * LATTICE_SIZE);
            let y = (idx / LATTICE_SIZE) % LATTICE_SIZE;
            let z = idx % LATTICE_SIZE;

            let spin = lattice.spins[[x, y, z]];
            let sz = spin.expectation_sz();

            // Map sz to a color
            let color = if sz > 0.0 {
                Point3::new(1.0 - sz as f32, 0.0, sz as f32) // Shades of blue
            } else {
                Point3::new(0.0, -sz as f32, 1.0 + sz as f32) // Shades of red
            };
            cube.set_color(color.x, color.y, color.z);
        }

        // Visualize the photon
        if let Some(photon) = &lattice.photon {
            let mut sphere = window.add_sphere(0.3);
            sphere.set_local_translation(
                Vector3::new(
                    photon.position[0] as f32,
                    photon.position[1] as f32,
                    photon.position[2] as f32,
                )
                .into(),
            );
            sphere.set_color(1.0, 1.0, 0.0); // Yellow color for the photon
        }
    }
}