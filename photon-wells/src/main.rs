use ndarray::prelude::*;
use num_complex::Complex;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;

/// Constants
const C: f64 = 299792458.0;          // Speed of light in vacuum (m/s)
const HBAR: f64 = 1.0545718e-34;     // Reduced Planck constant (J·s)
const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
const MU_0: f64 = 1.2566370614e-6;   // Vacuum permeability (H/m)
const G_A_GAMMA: f64 = 1e-13;        // Axion-photon coupling constant (1/GeV)
const AXION_MASS: f64 = 1e-5;        // Axion mass (eV)
const LATTICE_SIZE: usize = 20;      // Size of the lattice (20x20x20)
const TIME_STEPS: usize = 10000;      // Number of time steps
const DELTA_T: f64 = 1e-7;          // Time step (s)
const REFRACTIVE_INDEX_BASE: f64 = 1.0;  // Base refractive index
const NOISE_LEVEL: f64 = 1e-5;       // Noise level for susceptibility

const TUNNELING_PROB: f64 = 1e-4;    // Probability of tunneling event
const EVENT_HORIZON_N: f64 = 1.5;     // Threshold refractive index

// Neutrino parameters
const NEUTRINO_HALF_LIFE: f64 = 1e-15; // Example half-life in seconds
                                      // Probability of decay each timestep: p = 1 - exp(-ln(2)*DELTA_T/NEUTRINO_HALF_LIFE)

// Movement probability
const NEUTRINO_MOVE_PROB: f64 = 0.1; // Probability neutrino moves to a neighbor each step

/// Structure representing the axion field
#[derive(Clone, Copy, Debug)]
struct AxionField {
    value: f64,         // Axion field value at a point
    gradient: [f64; 3], // Spatial gradient of the axion field
}

impl AxionField {
    fn refractive_index_modification(&self) -> f64 {
        G_A_GAMMA * self.value
    }
}

/// Structure representing a photonic state at a point in the lattice
#[derive(Clone, Copy, Debug)]
struct PhotonicState {
    electric_field: Complex<f64>, // Electric field amplitude
    magnetic_field: Complex<f64>, // Magnetic field amplitude
    phase: f64,                   // Phase of the photon
}

impl PhotonicState {
    fn new_with_amplitude(rng: &mut StdRng, amplitude: f64) -> Self {
        let phase = rng.gen_range(0.0..(2.0 * PI));
        PhotonicState {
            electric_field: Complex::new(amplitude * phase.cos(), amplitude * phase.sin()),
            magnetic_field: Complex::new(amplitude * phase.cos(), amplitude * phase.sin()),
            phase,
        }
    }

    fn update_state(&mut self, delta_phase: f64) {
        let new_phase = (self.phase + delta_phase) % (2.0 * PI);
        let amplitude = self.electric_field.norm();
        self.electric_field = Complex::new(amplitude * new_phase.cos(), amplitude * new_phase.sin());
        self.magnetic_field = self.electric_field;
        self.phase = new_phase;
    }
}

/// Representing a neutrino spinor state
#[derive(Clone, Copy, Debug)]
struct NeutrinoState {
    spinor: [Complex<f64>; 4],
    flavor: usize, // 0, 1, or 2 representing neutrino flavor
    lifetime: f64, // Time accumulated to track decay
    active: bool,  // If false, no neutrino present
}

impl NeutrinoState {
    fn new_random(rng: &mut StdRng) -> Self {
        let spinor = [
            Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
            Complex::new(rng.gen::<f64>(), rng.gen::<f64>()),
        ];
        let norm = (spinor.iter().map(|c| c.norm_sqr()).sum::<f64>()).sqrt();
        let spinor = spinor.map(|c| c / norm);
        let flavor = rng.gen_range(0..3);
        NeutrinoState {
            spinor,
            flavor,
            lifetime: 0.0,
            active: true,
        }
    }

    fn update(
        &mut self,
        neighbors: &[NeutrinoState],
        axion: &AxionField,
        photon_intensity: f64,
        rng: &mut StdRng,
    ) {
        if !self.active {
            return;
        }
        let rotation_angle = (axion.value * photon_intensity).sin();
        let flavor = self.flavor;
        let new_flavor = ((flavor as f64 + rotation_angle) as usize) % 3;
        self.flavor = new_flavor;

        let flavor_count = neighbors.iter().filter(|n| n.active && n.flavor == self.flavor).count();
        if flavor_count > neighbors.len() / 2 {
            self.spinor[3] = -self.spinor[3];
        }

        let norm = (self.spinor.iter().map(|c| c.norm_sqr()).sum::<f64>()).sqrt();
        if norm != 0.0 {
            for s in &mut self.spinor {
                *s = *s / norm;
            }
        }
    }

    fn entangle_with_photon(&mut self, photon: &mut PhotonicState) {
        if !self.active {
            return;
        }
        let photon_amp = photon.electric_field.norm();
        let factor = Complex::new(photon_amp.cos(), photon_amp.sin());
        for s in &mut self.spinor {
            *s = (*s + factor) / Complex::new(2.0_f64.sqrt(), 0.0);
        }

        let new_amp = photon_amp / 2.0_f64.sqrt();
        photon.electric_field = Complex::new(new_amp * photon.phase.cos(), new_amp * photon.phase.sin());
        photon.magnetic_field = photon.electric_field;
    }

    fn decay_check(&mut self, rng: &mut StdRng) -> bool {
        // Probability of decay in this timestep:
        // p = 1 - exp(-ln(2)*DELTA_T/NEUTRINO_HALF_LIFE)
        let decay_prob = 1.0 - (- (std::f64::consts::LN_2 * DELTA_T / NEUTRINO_HALF_LIFE)).exp();
        rng.gen::<f64>() < decay_prob
    }

    fn decay_into_daughters(&mut self, rng: &mut StdRng) -> Vec<NeutrinoState> {
        // For simplicity, let's say one daughter neutrino is produced.
        // More complex models might produce multiple with different flavors.
        let mut daughter = NeutrinoState::new_random(rng);
        daughter.flavor = (self.flavor + 1) % 3; // shift flavor
        // Lower amplitude daughter: scale spinor slightly
        for s in &mut daughter.spinor {
            *s = *s * 0.5;
        }
        vec![daughter]
    }
}

/// Photonic + Neutrino Lattice
struct SimulationLattice {
    size: usize,
    photons: Array3<PhotonicState>,
    axion_field: Array3<AxionField>,
    neutrinos: Array3<NeutrinoState>,
}

impl SimulationLattice {
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(42);

        let center = (size / 2) as f64;
        let axion_field = Array3::from_shape_fn((size, size, size), |(x, y, z)| {
            let position = [x as f64 - center, y as f64 - center, z as f64 - center];
            let distance = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt();
            let value = (AXION_MASS * distance / center).sin();
            let gradient = [
                AXION_MASS * position[0].cos() / center,
                AXION_MASS * position[1].cos() / center,
                AXION_MASS * position[2].cos() / center,
            ];
            AxionField { value, gradient }
        });

        let radius = (size as f64) / 4.0;
        let photons = Array3::from_shape_fn((size, size, size), |(x, y, z)| {
            let position = [x as f64 - center, y as f64 - center, z as f64 - center];
            let dist = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt();
            let amplitude = if dist < radius { 1e-10 } else { 1.0 };
            PhotonicState::new_with_amplitude(&mut rng, amplitude)
        });

        let neutrinos = Array3::from_shape_fn((size, size, size), |_| {
            // Initially, let's say all sites have a neutrino
            let mut nu = NeutrinoState::new_random(&mut rng);
            // If desired, one could choose to have some sites empty:
            // For now, let's keep them all active
            nu
        });

        SimulationLattice {
            size,
            photons,
            axion_field,
            neutrinos,
        }
    }

    fn get_neighbors_neutrino(&self, x: usize, y: usize, z: usize) -> Vec<NeutrinoState> {
        let mut neighbors = Vec::new();
        let directions: [(isize, isize, isize); 6] = [
            (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)
        ];
        for (dx,dy,dz) in directions {
            let nx = x as isize + dx;
            let ny = y as isize + dy;
            let nz = z as isize + dz;
            if nx >= 0 && (nx as usize) < self.size &&
               ny >= 0 && (ny as usize) < self.size &&
               nz >= 0 && (nz as usize) < self.size {
                neighbors.push(self.neutrinos[[nx as usize, ny as usize, nz as usize]]);
            }
        }
        neighbors
    }

    fn evolve(&mut self) {
        let mut rng = StdRng::seed_from_u64(24);

        for _ in 0..TIME_STEPS {
            let photons_copy = self.photons.clone();
            let neutrinos_copy = self.neutrinos.clone();

            // We'll store changes in neutrinos here, including movements and decays
            let mut new_neutrinos = neutrinos_copy.clone();

            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        let photon = photons_copy[[x, y, z]];
                        let axion = self.axion_field[[x, y, z]];
                        let mut neutrino = neutrinos_copy[[x, y, z]];

                        if neutrino.active {
                            let delta_n = axion.refractive_index_modification();
                            let refractive_index = REFRACTIVE_INDEX_BASE + delta_n;

                            // Photon phase update
                            let omega = 2.0 * PI * C / (DELTA_T * refractive_index);
                            let delta_phase = omega * DELTA_T;
                            let mut new_photon = photon;
                            new_photon.update_state(delta_phase);

                            // Add noise
                            let noise_phase = rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
                            new_photon.update_state(noise_phase);

                            // Horizon-like emission
                            if refractive_index > EVENT_HORIZON_N {
                                new_photon = PhotonicState::new_with_amplitude(&mut rng, 1.0);
                            }

                            // Update neutrino (oscillation)
                            let neighbors = self.get_neighbors_neutrino(x, y, z);
                            let photon_intensity = new_photon.electric_field.norm_sqr();
                            neutrino.update(&neighbors, &axion, photon_intensity, &mut rng);

                            // Tunneling entanglement
                            if rng.gen::<f64>() < TUNNELING_PROB {
                                neutrino.entangle_with_photon(&mut new_photon);
                            }

                            // Decay check
                            neutrino.lifetime += DELTA_T;
                            if neutrino.decay_check(&mut rng) {
                                // Decay: remove current neutrino, produce daughters
                                neutrino.active = false;
                                // Place daughter neutrinos in neighboring cells if possible
                                let directions: [(isize, isize, isize); 6] = [
                                    (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)
                                ];
                                let daughters = neutrino.decay_into_daughters(&mut rng);
                                for daughter in daughters {
                                    // Attempt to place daughter in a random neighbor
                                    let (dx,dy,dz) = directions[rng.gen_range(0..6)];
                                    let nx = x as isize + dx;
                                    let ny = y as isize + dy;
                                    let nz = z as isize + dz;
                                    if nx >= 0 && (nx as usize) < self.size &&
                                       ny >= 0 && (ny as usize) < self.size &&
                                       nz >= 0 && (nz as usize) < self.size {
                                        // If target cell is empty or we allow multiple neutrinos
                                        // Here we overwrite for simplicity
                                        new_neutrinos[[nx as usize, ny as usize, nz as usize]] = daughter;
                                    }
                                }
                            } else {
                                // Neutrino survives. Consider movement.
                                if rng.gen::<f64>() < NEUTRINO_MOVE_PROB {
                                    let directions: [(isize, isize, isize); 6] = [
                                        (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)
                                    ];
                                    let (dx,dy,dz) = directions[rng.gen_range(0..6)];
                                    let nx = x as isize + dx;
                                    let ny = y as isize + dy;
                                    let nz = z as isize + dz;
                                    if nx >= 0 && (nx as usize) < self.size &&
                                       ny >= 0 && (ny as usize) < self.size &&
                                       nz >= 0 && (nz as usize) < self.size {
                                        // Move neutrino
                                        new_neutrinos[[x,y,z]].active = false;
                                        new_neutrinos[[nx as usize, ny as usize, nz as usize]] = neutrino;
                                    } else {
                                        // Can't move, stays in place
                                        new_neutrinos[[x,y,z]] = neutrino;
                                    }
                                } else {
                                    // No movement, just stay
                                    new_neutrinos[[x,y,z]] = neutrino;
                                }
                            }

                            self.photons[[x, y, z]] = new_photon;
                        } else {
                            // No neutrino active, just copy photon
                            self.photons[[x, y, z]] = photon;
                        }
                    }
                }
            }

            self.neutrinos = new_neutrinos;
        }
    }

    fn calculate_average_intensity(&self) {
        let mut total_intensity = 0.0;
        for photon in self.photons.iter() {
            total_intensity += photon.electric_field.norm_sqr();
        }
        let avg_intensity = total_intensity / (self.size.pow(3) as f64);
        println!("Average Electric Field Intensity: {}", avg_intensity);
    }

    fn analyze_noise(&self) {
        let intensities: Vec<f64> = self.photons.iter().map(|p| p.electric_field.norm_sqr()).collect();
        let mean_intensity = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let variance = intensities
            .iter()
            .map(|i| (i - mean_intensity).powi(2))
            .sum::<f64>() / intensities.len() as f64;
        println!("Intensity Variance (Noise Susceptibility): {}", variance);
    }

    fn visualize_slice(&self, filename: &str) {
        use plotters::prelude::*;
        let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let max_intensity = self
            .photons
            .iter()
            .map(|p| p.electric_field.norm_sqr())
            .fold(f64::NAN, f64::max);

        let z = self.size / 2;

        let mut chart = ChartBuilder::on(&root)
            .caption("Electric Field Intensity Slice", ("sans-serif", 20))
            .build_cartesian_2d(0..self.size, 0..self.size)
            .unwrap();
        chart.configure_mesh().draw().unwrap();

        for x in 0..self.size {
            for y in 0..self.size {
                let photon = self.photons[[x, y, z]];
                let intensity = photon.electric_field.norm_sqr();
                let color_value = if max_intensity.is_normal() {
                    (intensity / max_intensity * 255.0) as u8
                } else {
                    0
                };
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

fn main() {
    let mut lattice = SimulationLattice::new(LATTICE_SIZE);
    lattice.calculate_average_intensity();
    lattice.analyze_noise();
    lattice.visualize_slice("initial_intensity.png");

    for x in 0..100 {
    lattice.evolve();
    }

    lattice.calculate_average_intensity();
    lattice.analyze_noise();
    lattice.visualize_slice("final_intensity.png");

    println!("Simulation complete with neutrino decay, free movement, and fallback von Neumann algebraic reasoning.");
}