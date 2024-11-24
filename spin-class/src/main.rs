use num_complex::Complex;
use ndarray::prelude::*;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::Normal;
use plotters::prelude::*;

/// Constants
const HBAR: f64 = 1.0545718e-34; // Reduced Planck constant in J·s
const MU_B: f64 = 9.274009994e-24; // Bohr magneton in J/T
const KB: f64 = 1.380649e-23; // Boltzmann constant in J/K
const J_EXCHANGE: f64 = 1e-21; // Exchange interaction energy in J
const LATTICE_SIZE: usize = 20; // Size of the lattice (10x10x10)
const TIME_STEPS: usize = 1000; // Number of time steps
const DELTA_T: f64 = 1e-22; // Time step in seconds
const TEMPERATURE: f64 = 100000000000.0 * (1.0/137.0); // Temperature in Kelvin
const EXTERNAL_FIELD: f64 = -1000000090000000000.0 * (137.0/1.0);// -100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009000000000000000000.0 * (1.0/137.0); // External magnetic field in Tesla

// Added CMB temperature
const CMB_TEMPERATURE: f64 = 2.725; // Cosmic Microwave Background temperature in Kelvin

/// Spinor struct representing a quantum spin state
#[derive(Clone, Copy, Debug)]
struct Spinor {
    up: Complex<f64>,
    down: Complex<f64>,
}

impl Spinor {
    /// Initialize a spinor pointing in a random direction (to reflect permutation symmetry)
    fn new_random(rng: &mut StdRng) -> Self {
        // Random angles for spherical coordinates
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let phi = rng.gen_range(0.0..(2.0 * std::f64::consts::PI));

        // Convert to spinor components
        let up = Complex::new(theta.cos() / 2.0, 0.0);
        let down = Complex::new(theta.sin() / 2.0 * phi.cos(), theta.sin() / 2.0 * phi.sin());

        Spinor { up, down }
    }

    /// Normalize the spinor
    fn normalize(&mut self) {
        let norm = (self.up.norm_sqr() + self.down.norm_sqr()).sqrt();
        self.up /= norm;
        self.down /= norm;
    }

    /// Expectation value of Sx
    fn expectation_sx(&self) -> f64 {
        let sx = self.up.conj() * self.down + self.down.conj() * self.up;
        sx.re
    }

    /// Expectation value of Sy
    fn expectation_sy(&self) -> f64 {
        let sy: Complex<f64> = Complex::<f64>::i() * (self.up.conj() * self.down - self.down.conj() * self.up);
        sy.re
    }

    /// Expectation value of Sz
    fn expectation_sz(&self) -> f64 {
        let sz = self.up.conj() * self.up - self.down.conj() * self.down;
        sz.re
    }
}

/// Lattice struct representing the 3D lattice of spins
struct Lattice {
    spins: Array3<Spinor>,
    size: usize,
}

impl Lattice {
    /// Initialize a new lattice with spins in random orientations to reflect permutation symmetry
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(0);
        let spins = Array3::from_shape_fn((size, size, size), |_| {
            let mut spin = Spinor::new_random(&mut rng);
            spin.normalize();
            spin
        });
        Lattice { spins, size }
    }

    /// Simulate the evolution of the lattice over time
    fn evolve(&mut self) {
        let mut rng = StdRng::seed_from_u64(0);
        for _ in 0..TIME_STEPS {
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
                        let thermal_std = (2.0 * KB * TEMPERATURE / (MU_B)).sqrt(); // Thermal field standard deviation
                        let normal_dist = Normal::new(0.0, thermal_std).unwrap();
                        let thermal_field = [
                            rng.sample(normal_dist),
                            rng.sample(normal_dist),
                            rng.sample(normal_dist),
                        ];

                        // CMB field fluctuations (in Tesla)
                        let cmb_field = self.calculate_cmb_field(&mut rng);

                        // External magnetic field (set to zero in this simulation)
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
                        let delta = -DELTA_T; // Forward time evolution
                        let exponent = [
                            [Complex::new(h11, 0.0), h12],
                            [h21, Complex::new(h22, 0.0)],
                        ];
                        let exponent = matrix_scalar_multiply(&exponent, Complex::new(0.0, -delta / HBAR));

                        // Exponentiate the Hamiltonian matrix using Padé approximant for better numerical stability
                        let u_matrix = matrix_exponential(&exponent);

                        // Apply the time evolution operator
                        let new_up = u_matrix[0][0] * spin.up + u_matrix[0][1] * spin.down;
                        let new_down = u_matrix[1][0] * spin.up + u_matrix[1][1] * spin.down;

                        let mut new_spin = Spinor {
                            up: new_up,
                            down: new_down,
                        };

                        new_spin.normalize();

                        self.spins[[x, y, z]] = new_spin;
                    }
                }
            }
        }
    }

    /// Calculate the CMB field fluctuations for a spin
    fn calculate_cmb_field(&self, rng: &mut StdRng) -> [f64; 3] {
        // The CMB photons interact weakly, but we can model their effect as an additional thermal field
        let cmb_std = (2.0 * KB * CMB_TEMPERATURE / (MU_B)).sqrt(); // Standard deviation due to CMB
        let normal_dist = Normal::new(0.0, cmb_std).unwrap();
        [
            rng.sample(normal_dist),
            rng.sample(normal_dist),
            rng.sample(normal_dist),
        ]
    }

    /// Get the neighboring spins for a given position
    fn get_neighbors(&self, x: usize, y: usize, z: usize) -> Vec<Spinor> {
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

    /// Calculate the average magnetization of the lattice
    fn calculate_magnetization(&self) -> f64 {
        let mut total_sz = 0.0;
        for spin in self.spins.iter() {
            total_sz += spin.expectation_sz();
        }
        total_sz / (self.size * self.size * self.size) as f64
    }

    /// Calculate the Heisenberg uncertainty spread
    fn calculate_uncertainty(&self) -> f64 {
        let mut delta_sx_squared = 0.0;
        let mut delta_sy_squared = 0.0;
        let mut avg_sx = 0.0;
        let mut avg_sy = 0.0;
        for spin in self.spins.iter() {
            avg_sx += spin.expectation_sx();
            avg_sy += spin.expectation_sy();
        }
        avg_sx /= (self.size * self.size * self.size) as f64;
        avg_sy /= (self.size * self.size * self.size) as f64;
        for spin in self.spins.iter() {
            delta_sx_squared += (spin.expectation_sx() - avg_sx).powi(2);
            delta_sy_squared += (spin.expectation_sy() - avg_sy).powi(2);
        }
        delta_sx_squared /= (self.size * self.size * self.size) as f64;
        delta_sy_squared /= (self.size * self.size * self.size) as f64;
        (delta_sx_squared.sqrt()) * (delta_sy_squared.sqrt())
    }

    /// Plot a 2D slice of the magnetization
    fn plot_magnetization_slice(&self, filename: &str) {
        let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let max_value = 1.0;
        let min_value = -1.0;

        let mut chart = ChartBuilder::on(&root)
            .caption("Magnetization Slice", ("sans-serif", 20))
            .build_cartesian_2d(0..self.size, 0..self.size)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        let z = self.size / 2; // Plot the middle slice

        for x in 0..self.size {
            for y in 0..self.size {
                let spin = self.spins[[x, y, z]];
                let sz = spin.expectation_sz();
                let color_value = ((sz - min_value) / (max_value - min_value) * 255.0) as u8;
                chart
                    .draw_series(std::iter::once(Rectangle::new(
                        [(x, y), (x + 1, y + 1)],
                        RGBColor(color_value, 0, 255 - color_value).filled(),
                    )))
                    .unwrap();
            }
        }
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
        let coeff_num = factorial(n) * factorial(n - k);
        let coeff_den = factorial(n + k) * factorial(n - k);

        let coeff_num = Complex::new(coeff_num as f64, 0.0);
        let coeff_den = Complex::new(coeff_den as f64, 0.0);

        let term_num = matrix_scalar_multiply(&a_power, coeff_num);
        let term_den = matrix_scalar_multiply(&a_power, coeff_den);

        numerator = matrix_add(&numerator, &term_num);
        denominator = matrix_add(&denominator, &term_den);
    }

    // Inverse of denominator
    let det = denominator[0][0] * denominator[1][1] - denominator[0][1] * denominator[1][0];
    let inv_det = det.inv();
    let inverse_denominator = [
        [
            denominator[1][1] * inv_det,
            -denominator[0][1] * inv_det,
        ],
        [
            -denominator[1][0] * inv_det,
            denominator[0][0] * inv_det,
        ],
    ];

    matrix_multiply(&inverse_denominator, &numerator)
}

fn factorial(n: usize) -> usize {
    (1..=n).product()
}

fn main() {
    // Initialize the lattice with random spins to reflect permutation symmetry
    let mut lattice = Lattice::new(LATTICE_SIZE);

    // Initial magnetization and uncertainty
    let initial_magnetization = lattice.calculate_magnetization();
    let initial_uncertainty = lattice.calculate_uncertainty();
    println!("Initial Magnetization: {}", initial_magnetization);
    println!("Initial Uncertainty Spread: {}", initial_uncertainty);

    // Plot initial magnetization
    lattice.plot_magnetization_slice("initial_magnetization.png");

    // Evolve the lattice forward in time
    lattice.evolve();
        // Evolve the lattice forward in time
        lattice.evolve();
            // Evolve the lattice forward in time
    lattice.evolve();
        // Evolve the lattice forward in time
        lattice.evolve();

    // Magnetization and uncertainty after evolution
    let final_magnetization = lattice.calculate_magnetization();
    let final_uncertainty = lattice.calculate_uncertainty();
    println!("Final Magnetization: {}", final_magnetization);
    println!("Final Uncertainty Spread: {}", final_uncertainty);

    // Plot magnetization after evolution
    lattice.plot_magnetization_slice("final_magnetization.png");
}
