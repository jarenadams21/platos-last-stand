// Cargo.toml dependencies:
// [dependencies]
// num-complex = "0.4"
// ndarray = "0.15"
// rand = "0.8"
// plotters = "0.3"
// plotters-backend = "0.3"
// rand_distr = "0.4"

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
const LATTICE_SIZE: usize = 20; // Size of the lattice (20x20x20)
const TIME_STEPS: usize = 1000; // Number of time steps
const DELTA_T: f64 = 1e-17; // Time step in seconds
const TEMPERATURE: f64 = 1000000000.0; // Temperature in Kelvin
const EXTERNAL_FIELD: f64 = 10000.0; // External magnetic field in Tesla

// Added CMB temperature
const CMB_TEMPERATURE: f64 = 2.725; // Cosmic Microwave Background temperature in Kelvin

/// Spinor struct representing a quantum spin state
#[derive(Clone, Copy, Debug)]
struct Spinor {
    up: Complex<f64>,
    down: Complex<f64>,
}

impl Spinor {
    /// Initialize a spinor pointing up
    fn new_up() -> Self {
        Spinor {
            up: Complex::new(1.0, 0.0),
            down: Complex::new(0.0, 0.0),
        }
    }

    /// Initialize a spinor pointing down
    fn new_down() -> Self {
        Spinor {
            up: Complex::new(0.0, 0.0),
            down: Complex::new(1.0, 0.0),
        }
    }

    /// Normalize the spinor
    fn normalize(&mut self) {
        let norm = (self.up.norm_sqr() + self.down.norm_sqr()).sqrt();
        self.up /= norm;
        self.down /= norm;
    }

    /// Expectation value of Sx
    fn expectation_sx(&self) -> f64 {
        let sx: Complex<f64> = self.up.conj() * self.down + self.down.conj() * self.up;
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

/// Energy struct to keep track of energy exchanges
#[derive(Clone, Copy, Debug)]
struct Energy {
    total_energy: f64,
}

impl Energy {
    fn new() -> Self {
        Energy { total_energy: 0.0 }
    }

    fn add_energy(&mut self, delta_e: f64) {
        self.total_energy += delta_e;
    }
}

/// Lattice struct representing the 3D lattice of spins
struct Lattice {
    spins: Array3<Spinor>,
    size: usize,
    energy: Energy,
}

impl Lattice {
    /// Initialize a new lattice with all spins pointing up
    fn new(size: usize) -> Self {
        let spin_up = Spinor::new_up();
        let spins = Array3::from_elem((size, size, size), spin_up);
        Lattice {
            spins,
            size,
            energy: Energy::new(),
        }
    }

    /// Apply an external magnetic field in the center region (observer effect)
    fn apply_external_field(&mut self) {
        let center = self.size / 2;
        let radius = self.size / 5; // Define the non-magnetic sphere radius
        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    let dx = x as isize - center as isize;
                    let dy = y as isize - center as isize;
                    let dz = z as isize - center as isize;
                    let distance = ((dx * dx + dy * dy + dz * dz) as f64).sqrt();
                    if distance < radius as f64 {
                        // Flip the spins in the center region
                        self.spins[[x, y, z]] = Spinor::new_down();
                    }
                }
            }
        }
    }

    /// Simulate the evolution of the lattice over time
    fn evolve(&mut self, forward: bool) {
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

                        // Direction of time evolution
                        let time_factor = if forward { -1.0 } else { 1.0 };

                        // Magnetic moment of an electron spin (Bohr magneton)
                        const MU_S: f64 = MU_B;

                        // Hamiltonian matrix elements (in Joules)
                        let h11 = -0.5 * MU_S * total_field[2];
                        let h12 = -0.5 * MU_S * (total_field[0] - Complex::<f64>::i() * total_field[1]);
                        let h21 = -0.5 * MU_S * (total_field[0] + Complex::<f64>::i() * total_field[1]);
                        let h22 = 0.5 * MU_S * total_field[2];

                        // Time evolution operator: U = exp(-i * H * Δt / ħ)
                        let delta = time_factor * DELTA_T;
                        let exponent = [
                            [Complex::new(h11, 0.0), h12],
                            [h21, Complex::new(h22, 0.0)],
                        ];
                        let exponent = matrix_scalar_multiply(&exponent, Complex::new(0.0, -delta / HBAR));

                        // Exponentiate the Hamiltonian matrix using Taylor expansion
                        let identity = [
                            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
                            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
                        ];

                        let mut u_matrix = identity;
                        let mut term = identity;
                        let mut factorial = 1.0;

                        let n_terms = 20; // Number of terms in the Taylor series

                        for n in 1..n_terms {
                            factorial *= n as f64;
                            term = matrix_multiply(&term, &exponent);
                            let coeff = Complex::new(1.0 / factorial, 0.0);
                            u_matrix = matrix_add(&u_matrix, &matrix_scalar_multiply(&term, coeff));
                        }

                        // Apply the time evolution operator
                        let new_up = u_matrix[0][0] * spin.up + u_matrix[0][1] * spin.down;
                        let new_down = u_matrix[1][0] * spin.up + u_matrix[1][1] * spin.down;

                        let mut new_spin = Spinor {
                            up: new_up,
                            down: new_down,
                        };

                        new_spin.normalize();

                        // Calculate energy exchange due to CMB interaction
                        let energy_exchange = self.calculate_energy_exchange(&spin, &new_spin, &total_field);
                        self.energy.add_energy(energy_exchange);

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

    /// Calculate the energy exchange due to spin transition
    fn calculate_energy_exchange(&self, old_spin: &Spinor, new_spin: &Spinor, total_field: &[f64; 3]) -> f64 {
        // Energy difference ΔE = -μ · (B_new - B_old)
        let delta_b = total_field;
        let mu_s = MU_B;

        let s_old = [
            old_spin.expectation_sx(),
            old_spin.expectation_sy(),
            old_spin.expectation_sz(),
        ];
        let s_new = [
            new_spin.expectation_sx(),
            new_spin.expectation_sy(),
            new_spin.expectation_sz(),
        ];
        let delta_s = [
            s_new[0] - s_old[0],
            s_new[1] - s_old[1],
            s_new[2] - s_old[2],
        ];

        let delta_e = -mu_s * (delta_s[0] * delta_b[0] + delta_s[1] * delta_b[1] + delta_s[2] * delta_b[2]);
        delta_e
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

fn main() {
    // Initialize the lattice
    let mut lattice = Lattice::new(LATTICE_SIZE);

    // Apply external magnetic field (observer effect)
    lattice.apply_external_field();

    // Initial magnetization and uncertainty
    let initial_magnetization = lattice.calculate_magnetization();
    let initial_uncertainty = lattice.calculate_uncertainty();
    println!("Initial Magnetization: {}", initial_magnetization);
    println!("Initial Uncertainty Spread: {}", initial_uncertainty);

    // Plot initial magnetization
    lattice.plot_magnetization_slice("initial_magnetization.png");

    // Evolve the lattice forward in time
    lattice.evolve(true);

    // Magnetization and uncertainty after forward evolution
    let forward_magnetization = lattice.calculate_magnetization();
    let forward_uncertainty = lattice.calculate_uncertainty();
    println!("Forward Magnetization: {}", forward_magnetization);
    println!("Forward Uncertainty Spread: {}", forward_uncertainty);

    // Total energy after forward evolution
    println!("Total Energy after Forward Evolution: {}", lattice.energy.total_energy);

    // Plot magnetization after forward evolution
    lattice.plot_magnetization_slice("forward_magnetization.png");

    // Evolve the lattice backward in time
    lattice.evolve(false);

    // Magnetization and uncertainty after backward evolution
    let backward_magnetization = lattice.calculate_magnetization();
    let backward_uncertainty = lattice.calculate_uncertainty();
    println!("Backward Magnetization: {}", backward_magnetization);
    println!("Backward Uncertainty Spread: {}", backward_uncertainty);

    // Total energy after backward evolution
    println!("Total Energy after Backward Evolution: {}", lattice.energy.total_energy);

    // Plot magnetization after backward evolution
    lattice.plot_magnetization_slice("backward_magnetization.png");
}
