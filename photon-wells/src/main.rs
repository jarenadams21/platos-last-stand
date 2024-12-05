use ndarray::prelude::*;
use num_complex::Complex;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;

/// Constants
const C: f64 = 299792458.0; // Speed of light in vacuum (m/s)
const HBAR: f64 = 1.0545718e-34; // Reduced Planck constant (J·s)
const EPSILON_0: f64 = 8.854187817e-12; // Vacuum permittivity (F/m)
const MU_0: f64 = 1.2566370614e-6; // Vacuum permeability (H/m)
const G_A_GAMMA: f64 = 1e-13; // Axion-photon coupling constant (1/GeV)
const AXION_MASS: f64 = 1e-5; // Axion mass (eV)
const LATTICE_SIZE: usize = 20; // Size of the lattice (20x20x20)
const TIME_STEPS: usize = 1000; // Number of time steps
const DELTA_T: f64 = 1e-18; // Time step (s)
const TEMPERATURE: f64 = 2.725; // Temperature in Kelvin (CMB temperature)
const REFRACTIVE_INDEX_BASE: f64 = 1.0; // Base refractive index
const NOISE_LEVEL: f64 = 1e-5; // Noise level for susceptibility

/// Structure representing the axion field
#[derive(Clone, Copy, Debug)]
struct AxionField {
    value: f64,         // Axion field value at a point
    gradient: [f64; 3], // Spatial gradient of the axion field
}

impl AxionField {
    /// Calculate the axion-induced refractive index modification
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
    /// Initialize a photonic state with a random phase
    fn new_random(rng: &mut StdRng) -> Self {
        let phase = rng.gen_range(0.0..(2.0 * PI));
        PhotonicState {
            electric_field: Complex::new(phase.cos(), phase.sin()),
            magnetic_field: Complex::new(phase.cos(), phase.sin()),
            phase,
        }
    }

    /// Update the photonic state using trigonometric identities
    fn update_state(&mut self, delta_phase: f64) {
        // Use angle addition formula to compute new phase
        let new_phase = (self.phase + delta_phase) % (2.0 * PI);

        // Update electric and magnetic fields using exact trigonometric identities
        self.electric_field = Complex::new(new_phase.cos(), new_phase.sin());
        self.magnetic_field = Complex::new(new_phase.cos(), new_phase.sin());
        self.phase = new_phase;
    }
}

/// Structure representing the 3D photonic lattice
struct PhotonicLattice {
    size: usize,
    photons: Array3<PhotonicState>,
    axion_field: Array3<AxionField>,
}

impl PhotonicLattice {
    /// Initialize a new photonic lattice with random photonic states and axion field
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(42); // Seed for reproducibility

        // Initialize photonic states with random phases
        let photons = Array3::from_shape_fn((size, size, size), |_| {
            PhotonicState::new_random(&mut rng)
        });

        // Initialize axion field with a spatial gradient
        let axion_field = Array3::from_shape_fn((size, size, size), |(x, y, z)| {
            let center = (size / 2) as f64;
            let position = [
                x as f64 - center,
                y as f64 - center,
                z as f64 - center,
            ];
            let distance = (position[0].powi(2) + position[1].powi(2) + position[2].powi(2)).sqrt();
            let value = (AXION_MASS * distance / center).sin(); // Spatial variation using sine function
            let gradient = [
                AXION_MASS * position[0].cos() / center,
                AXION_MASS * position[1].cos() / center,
                AXION_MASS * position[2].cos() / center,
            ];
            AxionField { value, gradient }
        });

        PhotonicLattice {
            size,
            photons,
            axion_field,
        }
    }

    /// Simulate the evolution of the photonic lattice over time
    fn evolve(&mut self) {
        let mut rng = StdRng::seed_from_u64(24); // Seed for noise generation

        for _ in 0..TIME_STEPS {
            let photons_copy = self.photons.clone();

            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        // Get the current photonic state and axion field at this point
                        let photon = photons_copy[[x, y, z]];
                        let axion = self.axion_field[[x, y, z]];

                        // Calculate refractive index modification due to axion field
                        let delta_n = axion.refractive_index_modification();

                        // Total refractive index at this point
                        let refractive_index = REFRACTIVE_INDEX_BASE + delta_n;

                        // Calculate phase shift using exact trigonometric relations
                        // Phase shift: Δφ = (ω / c) * n * Δt
                        let omega = 2.0 * PI * C / (DELTA_T * refractive_index);
                        let delta_phase = omega * DELTA_T;

                        // Update the photonic state using trigonometric identities
                        let mut new_photon = photon;
                        new_photon.update_state(delta_phase);

                        // Noise due to thermal fluctuations and material imperfections
                        let noise_phase = rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
                        new_photon.update_state(noise_phase);

                        // Simulate Hawking radiation analogue by introducing an artificial event horizon
                        if refractive_index > 1.5 {
                            // At this threshold, photons are 'emitted' from the lattice
                            // Reset the photonic state to a new random state
                            self.photons[[x, y, z]] = PhotonicState::new_random(&mut rng);
                        } else {
                            self.photons[[x, y, z]] = new_photon;
                        }
                    }
                }
            }
        }
    }

    /// Calculate and print the average electric field intensity in the lattice
    fn calculate_average_intensity(&self) {
        let mut total_intensity = 0.0;
        for photon in self.photons.iter() {
            total_intensity += photon.electric_field.norm_sqr();
        }
        let average_intensity = total_intensity / (self.size * self.size * self.size) as f64;
        println!("Average Electric Field Intensity: {}", average_intensity);
    }

    /// Visualize a 2D slice of the lattice (e.g., electric field intensity)
    fn visualize_slice(&self, filename: &str) {
        use plotters::prelude::*;

        let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let max_intensity = self
            .photons
            .iter()
            .map(|p| p.electric_field.norm_sqr())
            .fold(0.0 / 0.0, f64::max); // Max intensity

        let z = self.size / 2; // Middle slice

        let mut chart = ChartBuilder::on(&root)
            .caption("Electric Field Intensity Slice", ("sans-serif", 20))
            .build_cartesian_2d(0..self.size, 0..self.size)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        for x in 0..self.size {
            for y in 0..self.size {
                let photon = self.photons[[x, y, z]];
                let intensity = photon.electric_field.norm_sqr();
                let color_value = (intensity / max_intensity * 255.0) as u8;
                chart
                    .draw_series(std::iter::once(Rectangle::new(
                        [(x, y), (x + 1, y + 1)],
                        RGBColor(color_value, 0, 255 - color_value).filled(),
                    )))
                    .unwrap();
            }
        }
    }

    /// Analyze noise susceptibility by measuring variance in field intensities
    fn analyze_noise(&self) {
        let intensities: Vec<f64> = self.photons.iter().map(|p| p.electric_field.norm_sqr()).collect();
        let mean_intensity = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let variance = intensities
            .iter()
            .map(|i| (i - mean_intensity).powi(2))
            .sum::<f64>()
            / intensities.len() as f64;
        println!("Intensity Variance (Noise Susceptibility): {}", variance);
    }
}

fn main() {
    // Initialize the photonic lattice
    let mut lattice = PhotonicLattice::new(LATTICE_SIZE);

    // Initial analysis
    lattice.calculate_average_intensity();
    lattice.analyze_noise();
    lattice.visualize_slice("initial_intensity.png");

    // Evolve the lattice over time
    lattice.evolve();

    // Post-evolution analysis
    lattice.calculate_average_intensity();
    lattice.analyze_noise();
    lattice.visualize_slice("final_intensity.png");

    println!("Simulation complete with enhanced algebraic relations.");
}