use rand::Rng;
use std::f64::consts::PI;

/// Constants
const HBAR: f64 = 1.0545718e-34; // Reduced Planck constant in J·s
const KB: f64 = 1.380649e-23;    // Boltzmann constant in J/K
const MU_B: f64 = 9.274009994e-24; // Bohr magneton in J/T
const J_EXCHANGE: f64 = 1e-21;   // Exchange energy constant in J
const LATTICE_SIZE: usize = 20;  // Size of the lattice (20x20x20)
const TIME_STEPS: usize = 100;   // Number of time steps
const TEMPERATURE: f64 = 300.0;  // Temperature in Kelvin
const EXTERNAL_FIELD: f64 = 1.0; // External magnetic field in Tesla

/// Spin struct representing a quantum spin state
#[derive(Clone, Copy)]
struct Spin {
    sx: f64,
    sy: f64,
    sz: f64,
}

impl Spin {
    /// Initialize a new spin pointing in the z-direction
    fn new_up() -> Self {
        Spin {
            sx: 0.0,
            sy: 0.0,
            sz: 1.0,
        }
    }

    /// Calculate the magnitude of the spin vector
    fn magnitude(&self) -> f64 {
        (self.sx.powi(2) + self.sy.powi(2) + self.sz.powi(2)).sqrt()
    }

    /// Normalize the spin vector
    fn normalize(&mut self) {
        let mag = self.magnitude();
        self.sx /= mag;
        self.sy /= mag;
        self.sz /= mag;
    }
}

/// Lattice struct representing the 3D lattice of spins
struct Lattice {
    spins: Vec<Vec<Vec<Spin>>>,
}

impl Lattice {
    /// Initialize a new lattice with all spins pointing up
    fn new() -> Self {
        let spin_up = Spin::new_up();
        let spins = vec![
            vec![vec![spin_up; LATTICE_SIZE]; LATTICE_SIZE];
            LATTICE_SIZE
        ];
        Lattice { spins }
    }

    /// Apply an external magnetic field in the center region
    fn apply_external_field(&mut self) {
        let center = LATTICE_SIZE / 2;
        let radius = LATTICE_SIZE / 5; // Define the non-magnetic sphere radius
        for x in 0..LATTICE_SIZE {
            for y in 0..LATTICE_SIZE {
                for z in 0..LATTICE_SIZE {
                    let dx = x as isize - center as isize;
                    let dy = y as isize - center as isize;
                    let dz = z as isize - center as isize;
                    let distance = ((dx * dx + dy * dy + dz * dz) as f64).sqrt();
                    if distance < radius as f64 {
                        // Flip the spins in the center region
                        self.spins[x][y][z] = Spin {
                            sx: 0.0,
                            sy: 0.0,
                            sz: -1.0,
                        };
                    }
                }
            }
        }
    }

    /// Simulate the evolution of the lattice over time
    fn evolve(&mut self) {
        let mut rng = rand::thread_rng();
        for _ in 0..TIME_STEPS {
            for x in 0..LATTICE_SIZE {
                for y in 0..LATTICE_SIZE {
                    for z in 0..LATTICE_SIZE {
                        let neighbors = self.get_neighbors(x, y, z);
                        let mut exchange_field = Spin { sx: 0.0, sy: 0.0, sz: 0.0 };
                        for neighbor in neighbors {
                            exchange_field.sx += neighbor.sx;
                            exchange_field.sy += neighbor.sy;
                            exchange_field.sz += neighbor.sz;
                        }
                        // Thermal fluctuations
                        let thermal_factor = (2.0 * rng.gen::<f64>() - 1.0)
                            * (2.0 * PI * KB * TEMPERATURE / HBAR);
                        // Effective field
                        let effective_field = Spin {
                            sx: J_EXCHANGE * exchange_field.sx + thermal_factor,
                            sy: J_EXCHANGE * exchange_field.sy + thermal_factor,
                            sz: J_EXCHANGE * exchange_field.sz
                                + thermal_factor
                                + MU_B * EXTERNAL_FIELD,
                        };
                        // Update spin using Landau-Lifshitz equation (simplified)
                        let current_spin = self.spins[x][y][z];
                        let cross_product = Spin {
                            sx: current_spin.sy * effective_field.sz
                                - current_spin.sz * effective_field.sy,
                            sy: current_spin.sz * effective_field.sx
                                - current_spin.sx * effective_field.sz,
                            sz: current_spin.sx * effective_field.sy
                                - current_spin.sy * effective_field.sx,
                        };
                        self.spins[x][y][z].sx += cross_product.sx * HBAR;
                        self.spins[x][y][z].sy += cross_product.sy * HBAR;
                        self.spins[x][y][z].sz += cross_product.sz * HBAR;
                        self.spins[x][y][z].normalize();
                    }
                }
            }
        }
    }

    /// Get the neighboring spins for a given position
    fn get_neighbors(&self, x: usize, y: usize, z: usize) -> Vec<Spin> {
        let mut neighbors = Vec::new();
        let positions = [
            (x.wrapping_sub(1), y, z),
            ((x + 1) % LATTICE_SIZE, y, z),
            (x, y.wrapping_sub(1), z),
            (x, (y + 1) % LATTICE_SIZE, z),
            (x, y, z.wrapping_sub(1)),
            (x, y, (z + 1) % LATTICE_SIZE),
        ];
        for &(nx, ny, nz) in &positions {
            if nx < LATTICE_SIZE && ny < LATTICE_SIZE && nz < LATTICE_SIZE {
                neighbors.push(self.spins[nx][ny][nz]);
            }
        }
        neighbors
    }

    /// Calculate the magnetization of the lattice
    fn calculate_magnetization(&self) -> f64 {
        let mut total_magnetization = 0.0;
        for x in 0..LATTICE_SIZE {
            for y in 0..LATTICE_SIZE {
                for z in 0..LATTICE_SIZE {
                    total_magnetization += self.spins[x][y][z].sz;
                }
            }
        }
        total_magnetization / (LATTICE_SIZE.pow(3) as f64)
    }

    /// Calculate the Heisenberg uncertainty spread
    fn calculate_uncertainty(&self) -> f64 {
        let mut delta_sx = 0.0;
        let mut delta_sy = 0.0;
        for x in 0..LATTICE_SIZE {
            for y in 0..LATTICE_SIZE {
                for z in 0..LATTICE_SIZE {
                    delta_sx += self.spins[x][y][z].sx.powi(2);
                    delta_sy += self.spins[x][y][z].sy.powi(2);
                }
            }
        }
        delta_sx = (delta_sx / (LATTICE_SIZE.pow(3) as f64)).sqrt();
        delta_sy = (delta_sy / (LATTICE_SIZE.pow(3) as f64)).sqrt();
        delta_sx * delta_sy
    }
}

fn main() {
    // Initialize the lattice
    let mut lattice = Lattice::new();

    // Apply external magnetic field (observer effect)
    lattice.apply_external_field();

    // Initial magnetization
    let initial_magnetization = lattice.calculate_magnetization();
    println!("Initial Magnetization: {}", initial_magnetization);

    // Initial uncertainty
    let initial_uncertainty = lattice.calculate_uncertainty();
    println!("Initial Uncertainty Spread: {}", initial_uncertainty);

    // Evolve the lattice over time
    lattice.evolve();

    // Final magnetization
    let final_magnetization = lattice.calculate_magnetization();
    println!("Final Magnetization: {}", final_magnetization);

    // Final uncertainty
    let final_uncertainty = lattice.calculate_uncertainty();
    println!("Final Uncertainty Spread: {}", final_uncertainty);
}
