use ndarray::prelude::*;
use num_complex::Complex;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;
use image::{RgbImage, Rgb};

/// Fundamental constants (SI units)
const C: f64 = 3.0e8;                 // Speed of light in m/s
const HBAR: f64 = 1.0545718e-34;      // Reduced Planck constant (J·s)
const G: f64 = 6.67430e-11;           // Gravitational constant (m^3 kg^-1 s^-2)
const K_B: f64 = 1.380649e-23;        // Boltzmann constant (J/K)

// Electromagnetic constants (not deeply used now, just placeholders)
const EPSILON_0: f64 = 8.854187817e-12; 
const MU_0: f64 = 1.2566370614e-6;   

// Field and mass parameters
const G_A_GAMMA: f64 = 1e-13;    // Axion-photon coupling constant (1/GeV), toy value
const AXION_MASS: f64 = 1e-5;    // Axion mass (eV), toy value

// Lattice parameters
const LATTICE_SIZE: usize = 14;
const TIME_STEPS: usize = 1000;
const DELTA_T: f64 = 1e-7;
const LATTICE_SPACING: f64 = 1e-3; // 1 mm spacing

// Black hole parameters
const RS: f64 = 1.0; // Schwarzschild radius in meters (example)
const REFRACTIVE_INDEX_BASE: f64 = 1.0;
const EVENT_HORIZON_N: f64 = 1.5;
const NOISE_LEVEL: f64 = 1e-5;
const TUNNELING_PROB: f64 = 1e-4;

// Neutrino parameters
const NEUTRINO_HALF_LIFE: f64 = 1e-15;
const NEUTRINO_MOVE_PROB: f64 = 0.1;

// Compute BH mass from RS
fn black_hole_mass(rs: f64) -> f64 {
    (rs * C.powi(2)) / (2.0 * G)
}

// Hawking temperature
fn hawking_temperature(m: f64) -> f64 {
    HBAR * C.powi(3) / (8.0 * PI * G * m * K_B)
}

// Structures
#[derive(Clone, Copy, Debug)]
struct AxionField {
    value: f64,
    gradient: [f64; 3],
}

impl AxionField {
    fn refractive_index_modification(&self) -> f64 {
        G_A_GAMMA * self.value
    }
}

#[derive(Clone, Copy, Debug)]
struct PhotonicState {
    electric_field: Complex<f64>,
    magnetic_field: Complex<f64>,
    phase: f64,
    position: [f64; 3],
}

impl PhotonicState {
    fn new_with_amplitude(rng: &mut StdRng, amplitude: f64, pos: [f64; 3]) -> Self {
        let phase = rng.gen_range(0.0..(2.0 * PI));
        PhotonicState {
            electric_field: Complex::new(amplitude * phase.cos(), amplitude * phase.sin() * -PI),
            magnetic_field: Complex::new(amplitude * phase.cos(), amplitude * phase.sin() * PI),
            phase,
            position: pos,
        }
    }

    fn update_state(&mut self, delta_phase: f64) {
        let new_phase = (self.phase + delta_phase) % (2.0 * PI);
        let amplitude = self.electric_field.norm();
        self.electric_field = Complex::new(amplitude * new_phase.cos(), amplitude * new_phase.sin());
        self.magnetic_field = self.electric_field;
        self.phase = new_phase;
    }

    fn intensity(&self) -> f64 {
        self.electric_field.norm_sqr()
    }
}

#[derive(Clone, Copy, Debug)]
struct NeutrinoState {
    spinor: [Complex<f64>; 4],
    flavor: usize,
    lifetime: f64,
    active: bool,
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
        self.flavor = ((self.flavor as f64 + rotation_angle) as usize) % 3;

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
        let decay_prob = 1.0 - (-(std::f64::consts::LN_2 * DELTA_T / NEUTRINO_HALF_LIFE)).exp();
        rng.gen::<f64>() < decay_prob
    }

    fn decay_into_daughters(&mut self, rng: &mut StdRng) -> Vec<NeutrinoState> {
        let mut daughter = NeutrinoState::new_random(rng);
        daughter.flavor = (self.flavor + 1) % 3;
        for s in &mut daughter.spinor {
            *s = *s * 0.5;
        }
        vec![daughter]
    }
}

struct SimulationLattice {
    size: usize,
    photons: Array3<PhotonicState>,
    axion_field: Array3<AxionField>,
    neutrinos: Array3<NeutrinoState>,
    center: f64,
    bh_mass: f64,
    bh_temperature: f64,
    rng: StdRng,
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
            let position = [x as f64, y as f64, z as f64];
            let dist = ((x as f64 - center).powi(2)
                      + (y as f64 - center).powi(2)
                      + (z as f64 - center).powi(2)).sqrt();
            let amplitude = if dist < radius { 1e-10 } else { 1.0 };
            PhotonicState::new_with_amplitude(&mut rng, amplitude, position)
        });

        let neutrinos = Array3::from_shape_fn((size, size, size), |_| {
            NeutrinoState::new_random(&mut rng)
        });

        let mass = black_hole_mass(RS);
        let temperature = hawking_temperature(mass);

        SimulationLattice {
            size,
            photons,
            axion_field,
            neutrinos,
            center,
            bh_mass: mass,
            bh_temperature: temperature,
            rng,
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

    fn effective_refractive_index(&self, x: usize, y: usize, z: usize) -> f64 {
        let position = [ (x as f64 - self.center)*LATTICE_SPACING,
                         (y as f64 - self.center)*LATTICE_SPACING,
                         (z as f64 - self.center)*LATTICE_SPACING ];
        let r = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt();
        let epsilon = 1e-6;
        REFRACTIVE_INDEX_BASE + RS/(r+epsilon)
    }

    // Sample photon amplitude from a thermal distribution at BH temperature.
    // For a rough approximation, let amplitude scale with sqrt(k_B T_H) to get a scale of fluctuation.
    // Actual photon energy distribution would be Planckian, but we simplify due to complexity.
    fn thermal_photon_amplitude(&mut self) -> f64 {
        let scale = (K_B * self.bh_temperature).sqrt();
        // Sample from a Maxwell-Boltzmann-like distribution:
        let u: f64 = self.rng.gen::<f64>();
        // Simple exponential decay distribution of amplitudes:
        (-u.ln()).sqrt() * scale
    }

    fn evolve(&mut self) {
        for _ in 0..TIME_STEPS {
            let photons_copy = self.photons.clone();
            let neutrinos_copy = self.neutrinos.clone();

            let mut new_neutrinos = neutrinos_copy.clone();

            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        let photon = photons_copy[[x, y, z]];
                        let axion = self.axion_field[[x, y, z]];
                        let mut neutrino = neutrinos_copy[[x, y, z]];

                        if neutrino.active {
                            let delta_n = axion.refractive_index_modification();
                            let refractive_index = self.effective_refractive_index(x, y, z) + delta_n;

                            // Update photon phase
                            let omega = 2.0 * PI * C / (DELTA_T * refractive_index);
                            let delta_phase = omega * DELTA_T;
                            let mut new_photon = photon;
                            new_photon.update_state(delta_phase);

                            // Add noise
                            let noise_phase = self.rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
                            new_photon.update_state(noise_phase);

                            // Horizon emission: sample amplitude from thermal distribution
                            if refractive_index > EVENT_HORIZON_N {
                                let amp = self.thermal_photon_amplitude();
                                new_photon = PhotonicState::new_with_amplitude(&mut self.rng, amp, photon.position);
                            }

                            // Neutrino update
                            let neighbors = self.get_neighbors_neutrino(x, y, z);
                            let photon_intensity = new_photon.intensity();
                            neutrino.update(&neighbors, &axion, photon_intensity, &mut self.rng);

                            // Tunneling entanglement
                            if self.rng.gen::<f64>() < TUNNELING_PROB {
                                neutrino.entangle_with_photon(&mut new_photon);
                            }

                            // Decay
                            neutrino.lifetime += DELTA_T;
                            if neutrino.decay_check(&mut self.rng) {
                                neutrino.active = false;
                                let daughters = neutrino.decay_into_daughters(&mut self.rng);
                                let directions: [(isize, isize, isize); 6] = [
                                    (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)
                                ];
                                for daughter in daughters {
                                    let (dx,dy,dz) = directions[self.rng.gen_range(0..6)];
                                    let nx = x as isize + dx;
                                    let ny = y as isize + dy;
                                    let nz = z as isize + dz;
                                    if nx >= 0 && (nx as usize) < self.size &&
                                       ny >= 0 && (ny as usize) < self.size &&
                                       nz >= 0 && (nz as usize) < self.size {
                                        new_neutrinos[[nx as usize, ny as usize, nz as usize]] = daughter;
                                    }
                                }
                            } else {
                                // Movement
                                if self.rng.gen::<f64>() < NEUTRINO_MOVE_PROB {
                                    let directions: [(isize, isize, isize); 6] = [
                                        (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)
                                    ];
                                    let (dx,dy,dz) = directions[self.rng.gen_range(0..6)];
                                    let nx = x as isize + dx;
                                    let ny = y as isize + dy;
                                    let nz = z as isize + dz;
                                    if nx >= 0 && (nx as usize) < self.size &&
                                       ny >= 0 && (ny as usize) < self.size &&
                                       nz >= 0 && (nz as usize) < self.size {
                                        new_neutrinos[[x,y,z]].active = false;
                                        new_neutrinos[[nx as usize, ny as usize, nz as usize]] = neutrino;
                                    } else {
                                        new_neutrinos[[x,y,z]] = neutrino;
                                    }
                                } else {
                                    new_neutrinos[[x,y,z]] = neutrino;
                                }
                            }

                            self.photons[[x, y, z]] = new_photon;
                        } else {
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
            total_intensity += photon.intensity();
        }
        let avg_intensity = total_intensity / (self.size.pow(3) as f64);
        println!("Average Electric Field Intensity: {}", avg_intensity);
    }

    fn analyze_noise(&self) {
        let intensities: Vec<f64> = self.photons.iter().map(|p| p.intensity()).collect();
        let mean_intensity = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let variance = intensities
            .iter()
            .map(|i| (i - mean_intensity).powi(2))
            .sum::<f64>() / intensities.len() as f64;
        println!("Intensity Variance (Noise Susceptibility): {}", variance);
    }

    fn sample_3d_vectors(&self) {
        let x = self.size / 2;
        let y = self.size / 2;
        let z = self.size / 2;
        let photon = self.photons[[x,y,z]];
        println!("Central Photon Electric Field Vector: ({:.3e}, {:.3e})",
                 photon.electric_field.re, photon.electric_field.im);
        println!("Central Photon Magnetic Field Vector: ({:.3e}, {:.3e})",
                 photon.magnetic_field.re, photon.magnetic_field.im);
        let axion = self.axion_field[[x,y,z]];
        println!("Axion Gradient at Center: {:?}", axion.gradient);
    }

    // Create a simple 2D visualization by projecting the photon intensity along one axis (e.g., z-axis).
    // We'll take the maximum intensity along z for each (x,y). This gives a simple 3D-to-2D projection.
    fn create_3d_visual(&self, filename: &str) {
        let mut image = RgbImage::new(self.size as u32, self.size as u32);

        // Max intensity projection along z
        for x in 0..self.size {
            for y in 0..self.size {
                let mut max_intensity = 0.0;
                for z in 0..self.size {
                    let intensity = self.photons[[x,y,z]].intensity();
                    if intensity > max_intensity {
                        max_intensity = intensity;
                    }
                }
                // Map intensity to color
                // Simple linear mapping: assume intensity ~0 to ~1 for demonstration.
                let val = (max_intensity * 255.0).min(255.0) as u8;
                image.put_pixel(x as u32, y as u32, Rgb([val, 0, 255 - val]));
            }
        }

        image.save(filename).unwrap();
        println!("3D visualization saved to {}", filename);
    }
}

fn main() {
    let mut lattice = SimulationLattice::new(LATTICE_SIZE);
    lattice.calculate_average_intensity();
    lattice.analyze_noise();
    lattice.sample_3d_vectors();

    lattice.evolve();

    println!("Evolutions by hand, done.");

    lattice.calculate_average_intensity();
    lattice.analyze_noise();

    lattice.create_3d_visual("3d_visualization.png");

    println!("Simulation complete with horizon-like behavior, thermal emission, and field interactions.");
}
