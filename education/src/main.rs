use ndarray::prelude::*;
use ndarray::azip;
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
const EPSILON_0: f64 = 8.854187817e-12;
const MU_0: f64 = 1.2566370614e-6;

// Field and mass parameters
const G_A_GAMMA: f64 = 1e-15;    // Axion-photon coupling constant
const AXION_MASS: f64 = 1e-10;    // Axion mass (eV), toy value

// Lattice parameters
const LATTICE_SIZE: usize = 20;
const TIME_STEPS: usize = 10;
const DELTA_T: f64 = 1e-22; // 15
const LATTICE_SPACING: f64 = 1e-6; // 1 mm spacing

// Black hole parameters
const RS: f64 = 1e-6; // Schwarzschild radius in meters
const REFRACTIVE_INDEX_BASE: f64 = 1.0;
const EVENT_HORIZON_N: f64 = 1.33; // (!) : Could be a sharp transition layer instead
const NOISE_LEVEL: f64 = 1e-4;
const TUNNELING_PROB: f64 = 1e-4;

// Neutrino parameters
const NEUTRINO_HALF_LIFE: f64 = 1e-10; // Unjustified
const NEUTRINO_MOVE_PROB: f64 = 0.2; // Unjustified

// Compute BH mass from RS
fn black_hole_mass(rs: f64) -> f64 {
    (rs * C.powi(2)) / (2.0 * G)
}

// Hawking temperature
fn hawking_temperature(m: f64) -> f64 {
    HBAR * C.powi(3) / (8.0 * PI * G * m * K_B)
}

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

/// Neutrino state
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
        if !self.active { return; }
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

    fn entangle_with_photon(&mut self, e_amp: &mut f64, phase: &mut f64) {
        if !self.active { return; }
        let photon_amp = *e_amp;
        let factor = Complex::new(photon_amp.cos(), photon_amp.sin());
        for s in &mut self.spinor {
            *s = (*s + factor) / Complex::new(2.0_f64.sqrt(), 0.0);
        }

        let new_amp = photon_amp / 2.0_f64.sqrt();
        *e_amp = new_amp;
        *phase = *phase % (2.0 * PI);
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

/// Simulation lattice with EM fields stored as separate arrays
struct SimulationLattice {
    size: usize,
    axion_field: Array3<AxionField>,
    neutrinos: Array3<NeutrinoState>,
    center: f64,
    bh_mass: f64,
    bh_temperature: f64,
    rng: StdRng,

    // Electric and magnetic fields (3D arrays)
    e_x: Array3<f64>,
    e_y: Array3<f64>,
    e_z: Array3<f64>,
    b_x: Array3<f64>,
    b_y: Array3<f64>,
    b_z: Array3<f64>,
}

impl SimulationLattice {
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(42);
        let center = (size / 2) as f64;
        let axion_field = Array3::from_shape_fn((size, size, size), |(x, y, z)| {
            let px = x as f64 - center;
            let py = y as f64 - center;
            let pz = z as f64 - center;
            let distance = (px*px + py*py + pz*pz).sqrt();
            let value = (AXION_MASS * distance).cos();
            let gradient = [
                -AXION_MASS * px.sin(),
                -AXION_MASS * py.sin(),
                -AXION_MASS * pz.sin() * PI * PI.sqrt(),
            ];
            AxionField { value, gradient }
        });

        let neutrinos = Array3::from_shape_fn((size, size, size), |_| {
            NeutrinoState::new_random(&mut rng)
        });

        let mass = black_hole_mass(RS);
        let temperature = hawking_temperature(mass);

        // Initialize E and B fields
        // Start with small random fields or uniform fields
        let e_x = Array3::from_shape_fn((size, size, size), |_| rng.gen_range(-1e-10..1e-10));
        let e_y = e_x.clone();
        let e_z = e_x.clone();
        let b_x = e_x.clone();
        let b_y = e_x.clone();
        let b_z = e_x.clone();

        SimulationLattice {
            size,
            axion_field,
            neutrinos,
            center,
            bh_mass: mass,
            bh_temperature: temperature,
            rng,
            e_x,
            e_y,
            e_z,
            b_x,
            b_y,
            b_z,
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
        let px = (x as f64 - self.center)*LATTICE_SPACING;
        let py = (y as f64 - self.center)*LATTICE_SPACING;
        let pz = (z as f64 - self.center)*LATTICE_SPACING;
        let r = (px*px + py*py + pz*pz).sqrt();
        let epsilon = 1e-6;
        REFRACTIVE_INDEX_BASE + RS/(r+epsilon)
    }

    fn thermal_photon_amplitude(&mut self) -> f64 {
        let scale = (K_B * self.bh_temperature).sqrt();
        let u: f64 = self.rng.gen::<f64>();
        (-u.ln()).sqrt() * scale
    }

    // Update E and B fields using a simple FDTD scheme
    // This is a rough approximation; boundary conditions and stability checks omitted for brevity
    fn update_em_fields(&mut self) {
        let dx = LATTICE_SPACING;
        let c = C;
        // Temporary arrays for updated fields
        let mut new_e_x = self.e_x.clone();
        let mut new_e_y = self.e_y.clone();
        let mut new_e_z = self.e_z.clone();
        let mut new_b_x = self.b_x.clone();
        let mut new_b_y = self.b_y.clone();
        let mut new_b_z = self.b_z.clone();

        // Curl computations (periodic boundaries for simplicity)
        for x in 0..self.size {
            let xp = (x+1) % self.size;
            let xm = (x+self.size-1) % self.size;
            for y in 0..self.size {
                let yp = (y+1) % self.size;
                let ym = (y+self.size-1) % self.size;
                for z in 0..self.size {
                    let zp = (z+1) % self.size;
                    let zm = (z+self.size-1) % self.size;

                    let refr_index = self.effective_refractive_index(x,y,z);
                    // Curl B
                    let curl_b_x = (self.b_z[[x,yp,z]] - self.b_z[[x,ym,z]])/(2.0*dx)
                                 - (self.b_y[[x,y,zp]] - self.b_y[[x,y,zm]])/(2.0*dx);

                    let curl_b_y = (self.b_x[[x,y,zp]] - self.b_x[[x,y,zm]])/(2.0*dx)
                                 - (self.b_z[[xp,y,z]] - self.b_z[[xm,y,z]])/(2.0*dx);

                    let curl_b_z = (self.b_y[[xp,y,z]] - self.b_y[[xm,y,z]])/(2.0*dx)
                                 - (self.b_x[[x,yp,z]] - self.b_x[[x,ym,z]])/(2.0*dx);

                    // Curl E
                    let curl_e_x = (self.e_z[[x,yp,z]] - self.e_z[[x,ym,z]])/(2.0*dx)
                                 - (self.e_y[[x,y,zp]] - self.e_y[[x,y,zm]])/(2.0*dx);

                    let curl_e_y = (self.e_x[[x,y,zp]] - self.e_x[[x,y,zm]])/(2.0*dx)
                                 - (self.e_z[[xp,y,z]] - self.e_z[[xm,y,z]])/(2.0*dx);

                    let curl_e_z = (self.e_y[[xp,y,z]] - self.e_y[[xm,y,z]])/(2.0*dx)
                                 - (self.e_x[[x,yp,z]] - self.e_x[[x,ym,z]])/(2.0*dx);

                    // Update E fields: dE/dt = c^2/n^2 curl B
                    let n_sq = refr_index*refr_index;
                    new_e_x[[x,y,z]] = self.e_x[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_x * PI;
                    new_e_y[[x,y,z]] = self.e_y[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_y * PI;
                    new_e_z[[x,y,z]] = self.e_z[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_z * PI;

                    // Update B fields: dB/dt = - curl E
                    new_b_x[[x,y,z]] = self.b_x[[x,y,z]] - DELTA_T * curl_e_x * -PI;
                    new_b_y[[x,y,z]] = self.b_y[[x,y,z]] - DELTA_T * curl_e_y * -PI;
                    new_b_z[[x,y,z]] = self.b_z[[x,y,z]] - DELTA_T * curl_e_z * -PI;
                }
            }
        }

        self.e_x = new_e_x;
        self.e_y = new_e_y;
        self.e_z = new_e_z;
        self.b_x = new_b_x;
        self.b_y = new_b_y;
        self.b_z = new_b_z;
    }

    fn evolve(&mut self) {
        for _ in 0..TIME_STEPS {
            let neutrinos_copy = self.neutrinos.clone();

            // First update EM fields
            self.update_em_fields();

            // Add noise and horizon emission
            for x in 0..self.size {
                for y in 0..self.size {
                    for z in 0..self.size {
                        let refr_index = self.effective_refractive_index(x,y,z) + self.axion_field[[x,y,z]].refractive_index_modification();
                        
                        // Add small noise to E fields
                        let noise = self.rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
                        self.e_x[[x,y,z]] += noise;
                        self.e_y[[x,y,z]] += noise;
                        self.e_z[[x,y,z]] += noise;

                        // Horizon emission if refr_index > EVENT_HORIZON_N
                        if refr_index > EVENT_HORIZON_N {
                            let amp = self.thermal_photon_amplitude();
                            // Random orientation for a photon-like state
                            let phase = self.rng.gen_range(0.0..2.0*PI);
                            self.e_x[[x,y,z]] = amp * phase.cos();
                            self.e_y[[x,y,z]] = amp * phase.sin();
                            self.e_z[[x,y,z]] = 0.0; // Simplify to 2D polarization
                        }

                        let mut neutrino = neutrinos_copy[[x,y,z]];
                        if neutrino.active {
                            let neighbors = self.get_neighbors_neutrino(x,y,z);
                            let intensity = (self.e_x[[x,y,z]].powi(2)
                                           + self.e_y[[x,y,z]].powi(2)
                                           + self.e_z[[x,y,z]].powi(2));
                            neutrino.update(&neighbors, &self.axion_field[[x,y,z]], intensity, &mut self.rng);

                            // Tunneling entanglement
                            if self.rng.gen::<f64>() < TUNNELING_PROB {
                                let mut phase = self.rng.gen_range(0.0..2.0*PI);
                                let mut e_amp = intensity.sqrt();
                                neutrino.entangle_with_photon(&mut e_amp, &mut phase);
                                // Update fields after entanglement
                                self.e_x[[x,y,z]] = e_amp * phase.cos();
                                self.e_y[[x,y,z]] = e_amp * phase.sin();
                                self.e_z[[x,y,z]] = 0.0;
                            }

                            // Decay check
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
                                        self.neutrinos[[nx as usize, ny as usize, nz as usize]] = daughter;
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
                                        self.neutrinos[[x,y,z]].active = false;
                                        self.neutrinos[[nx as usize, ny as usize, nz as usize]] = neutrino;
                                    } else {
                                        self.neutrinos[[x,y,z]] = neutrino;
                                    }
                                } else {
                                    self.neutrinos[[x,y,z]] = neutrino;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fn calculate_average_intensity(&self) {
        let mut total_intensity = 0.0;
        azip!((ex in &self.e_x, ey in &self.e_y, ez in &self.e_z) {
            total_intensity += ex*ex + ey*ey + ez*ez;
        });
        let avg_intensity = total_intensity / (self.size.pow(3) as f64);
        println!("Average Electric Field Intensity: {}", avg_intensity);
    }

    fn analyze_noise(&self) {
        let mut intensities = Vec::with_capacity(self.size.pow(3));
        azip!((ex in &self.e_x, ey in &self.e_y, ez in &self.e_z) {
            intensities.push(ex*ex + ey*ey + ez*ez);
        });

        let mean_intensity: f64 = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let variance = intensities.iter().map(|i| (i - mean_intensity).powi(2)).sum::<f64>() / intensities.len() as f64;
        println!("Intensity Variance (Noise Susceptibility): {}", variance);
    }

    fn sample_3d_vectors(&self) {
        let x = self.size / 2;
        let y = self.size / 2;
        let z = self.size / 2;
        println!("Central Electric Field Vector: ({:.3e}, {:.3e}, {:.3e})",
                 self.e_x[[x,y,z]], self.e_y[[x,y,z]], self.e_z[[x,y,z]]);
        println!("Central Magnetic Field Vector: ({:.3e}, {:.3e}, {:.3e})",
                 self.b_x[[x,y,z]], self.b_y[[x,y,z]], self.b_z[[x,y,z]]);
        let axion = self.axion_field[[x,y,z]];
        println!("Axion Gradient at Center: {:?}", axion.gradient);
    }

    fn create_3d_visual(&self, filename: &str) {
        let mut image = RgbImage::new(self.size as u32, self.size as u32);

        for x in 0..self.size {
            for y in 0..self.size {
                let mut max_intensity = 0.0;
                for z in 0..self.size {
                    let intensity = self.e_x[[x,y,z]].powi(2)
                                 + self.e_y[[x,y,z]].powi(2)
                                 + self.e_z[[x,y,z]].powi(2);
                    if intensity > max_intensity {
                        max_intensity = intensity;
                    }
                }
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

    for z in 0..1000 {
    lattice.evolve();
    lattice.calculate_average_intensity();
    }

    lattice.calculate_average_intensity();
    lattice.analyze_noise();

    lattice.create_3d_visual("3d_visualization.png");

    println!("Simulation complete with horizon-like behavior, Maxwellian EM fields, thermal emission, and field interactions. Researchers can adjust parameters for novel investigations.");
}
