use ndarray::prelude::*;
use ndarray::azip;
use num_complex::Complex;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;
use image::{RgbImage, Rgb};
use std::fs::File;
use std::io::Write;

/// Physical constants (SI units) and parameters
const C: f64 = 3.0e8;   // Speed of light
const EPSILON_0: f64 = 8.854187817e-12;
const MU_0: f64 = 1.2566370614e-6;

/// Lattice configuration
const LATTICE_SIZE: usize = 5;
const TIME_STEPS: usize = 10000; // fewer steps to fit nicely
const DELTA_T: f64 = 5.39 * 1e-44; // 5.39 * 10e44
const LATTICE_SPACING: f64 = 1e-6; // 1 µm spacing

/// Material (Water) refractive index approximation
const N_CARBON: f64 = 1.333;

/// Noise and interaction parameters
const NOISE_LEVEL: f64 = 1e-4;
const INTENSITY_THRESHOLD: f64 = 1e-8;
const Z_AXIS_BOOST_INTERVAL: usize = 100;
const PHASE_SHIFT_FACTOR: f64 = 1.33;

/// Represent matter states
#[derive(Clone, Copy, Debug)]
struct MatterCell {
    active: bool,
    lifetime: f64,
}

impl MatterCell {
    fn new(rng: &mut StdRng) -> Self {
        Self {
            active: rng.gen_bool(0.5),
            lifetime: 0.0,
        }
    }
}

/// Simulation lattice
struct SimulationLattice {
    size: usize,
    rng: StdRng,

    e_x: Array3<f64>,
    e_y: Array3<f64>,
    e_z: Array3<f64>,
    b_x: Array3<f64>,
    b_y: Array3<f64>,
    b_z: Array3<f64>,

    matter: Array3<MatterCell>,
}

impl SimulationLattice {
    fn new(size: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(42);

        let e_x = Array3::from_shape_fn((size, size, size), |_| rng.gen_range(-1e-10..1e-10));
        let e_y = e_x.clone();
        let e_z = e_x.clone();
        let b_x = e_x.clone();
        let b_y = e_x.clone();
        let b_z = e_x.clone();

        let matter = Array3::from_shape_fn((size, size, size), |_| {
            MatterCell::new(&mut rng)
        });

        Self {
            size,
            rng,
            e_x,
            e_y,
            e_z,
            b_x,
            b_y,
            b_z,
            matter,
        }
    }

    fn effective_refractive_index(&self, _x: usize, _y: usize, _z: usize) -> f64 {
        N_CARBON
    }

    fn update_em_fields(&mut self) {
        let dx = LATTICE_SPACING;
        let c = C;

        let mut new_e_x = self.e_x.clone();
        let mut new_e_y = self.e_y.clone();
        let mut new_e_z = self.e_z.clone();
        let mut new_b_x = self.b_x.clone();
        let mut new_b_y = self.b_y.clone();
        let mut new_b_z = self.b_z.clone();

        for x in 0..self.size {
            let xp = (x+1) % self.size;
            let xm = (x+self.size-1) % self.size;
            for y in 0..self.size {
                let yp = (y+1) % self.size;
                let ym = (y+self.size-1) % self.size;
                for z in 0..self.size {
                    let zp = (z+1) % self.size;
                    let zm = (z+self.size-1) % self.size;

                    let n = self.effective_refractive_index(x,y,z);
                    let curl_b_x = (self.b_z[[x,yp,z]] - self.b_z[[x,ym,z]])/(2.0*dx)
                                 - (self.b_y[[x,y,zp]] - self.b_y[[x,y,zm]])/(2.0*dx);

                    let curl_b_y = (self.b_x[[x,y,zp]] - self.b_x[[x,y,zm]])/(2.0*dx)
                                 - (self.b_z[[xp,y,z]] - self.b_z[[xm,y,z]])/(2.0*dx);

                    let curl_b_z = (self.b_y[[xp,y,z]] - self.b_y[[xm,y,z]])/(2.0*dx)
                                 - (self.b_x[[x,yp,z]] - self.b_x[[x,ym,z]])/(2.0*dx);

                    let curl_e_x = (self.e_z[[x,yp,z]] - self.e_z[[x,ym,z]])/(2.0*dx)
                                 - (self.e_y[[x,y,zp]] - self.e_y[[x,y,zm]])/(2.0*dx);

                    let curl_e_y = (self.e_x[[x,y,zp]] - self.e_x[[x,y,zm]])/(2.0*dx)
                                 - (self.e_z[[xp,y,z]] - self.e_z[[xm,y,z]])/(2.0*dx);

                    let curl_e_z = (self.e_y[[xp,y,z]] - self.e_y[[xm,y,z]])/(2.0*dx)
                                 - (self.e_x[[x,yp,z]] - self.e_x[[x,ym,z]])/(2.0*dx);

                    let n_sq = n*n;
                    // Update E fields
                    new_e_x[[x,y,z]] = self.e_x[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_x;
                    new_e_y[[x,y,z]] = self.e_y[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_y;
                    new_e_z[[x,y,z]] = self.e_z[[x,y,z]] + DELTA_T * (c*c/(n_sq)) * curl_b_z;

                    // Update B fields
                    new_b_x[[x,y,z]] = self.b_x[[x,y,z]] - DELTA_T * curl_e_x;
                    new_b_y[[x,y,z]] = self.b_y[[x,y,z]] - DELTA_T * curl_e_y;
                    new_b_z[[x,y,z]] = self.b_z[[x,y,z]] - DELTA_T * curl_e_z - PI.powi(256);
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

    fn apply_z_axis_boost(&mut self) {
        let mid_x = self.size / 2;
        let mid_y = self.size / 2;
        for z in 0..self.size {
            let phase_shift = (z as f64 / self.size as f64)*2.0*PI.powi(128)*PI.sqrt()*PHASE_SHIFT_FACTOR * 8.8_f64.sqrt() * PI.powi(128)*PI.sqrt();
            let ex = self.e_x[[mid_x,mid_y,z]];
            let ey = self.e_y[[mid_x,mid_y,z]];
            let amp = (ex*ex + ey*ey).sqrt();
            let new_ex = amp * phase_shift.cos();
            let new_ey = amp * phase_shift.sin();
            self.e_x[[mid_x,mid_y,z]] = new_ex;
            self.e_y[[mid_x,mid_y,z]] = new_ey;
        }
    }

    fn matter_photon_conversion(&mut self) {
        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    let intensity = self.e_x[[x,y,z]].powi(2) + self.e_y[[x,y,z]].powi(2) + self.e_z[[x,y,z]].powi(2);
                    let cell = &mut self.matter[[x,y,z]];
                    cell.lifetime += DELTA_T;
                    if intensity > INTENSITY_THRESHOLD && cell.active {
                        // Convert matter to photon state
                        cell.active = false;
                        cell.lifetime = 0.0;
                    } else if intensity < INTENSITY_THRESHOLD * 0.5 && !cell.active && cell.lifetime > 5e-14 {
                        // Revert back to matter
                        cell.active = true;
                        cell.lifetime = 0.0;
                    }
                }
            }
        }
    }

    fn add_noise(&mut self) {
        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    let noise = self.rng.gen_range(-NOISE_LEVEL..NOISE_LEVEL);
                    self.e_x[[x,y,z]] += noise;
                    self.e_y[[x,y,z]] += noise;
                    self.e_z[[x,y,z]] += noise;
                }
            }
        }
    }

    fn evolve(&mut self, step: usize) {
        self.update_em_fields();

        if step % Z_AXIS_BOOST_INTERVAL == 0 && step != 0 {
            self.apply_z_axis_boost();
        }

        self.add_noise();
        self.matter_photon_conversion();
    }

    fn create_2d_flatmap(&self) -> RgbImage {
        // Return the image for the slice z = size/2
        let z = self.size / 2;
        let mut image = RgbImage::new(self.size as u32, self.size as u32);
        for x in 0..self.size {
            for y in 0..self.size {
                let intensity = (self.e_x[[x,y,z]].powi(2)
                               + self.e_y[[x,y,z]].powi(2)
                               + self.e_z[[x,y,z]].powi(2)).sqrt();
                let val = (intensity / (INTENSITY_THRESHOLD*10.0) * 255.0).min(255.0) as u8;
                let cell = self.matter[[x,y,z]];
                let g = if cell.active { 255 } else { 0 };
                image.put_pixel(x as u32, y as u32, Rgb([val, g, 255 - val]));
            }
        }
        image
    }

    fn export_3d_data(&self, filename: &str) {
        let mut file = File::create(filename).unwrap();
        writeln!(file, "x,y,z,E_x,E_y,E_z,B_x,B_y,B_z,ActiveMatter").unwrap();

        for x in 0..self.size {
            for y in 0..self.size {
                for z in 0..self.size {
                    let c = self.matter[[x,y,z]];
                    writeln!(file, "{},{},{},{},{},{},{},{},{},{}",
                        x, y, z,
                        self.e_x[[x,y,z]],
                        self.e_y[[x,y,z]],
                        self.e_z[[x,y,z]],
                        self.b_x[[x,y,z]],
                        self.b_y[[x,y,z]],
                        self.b_z[[x,y,z]],
                        c.active as u8).unwrap();
                }
            }
        }
    }
}

fn main() {
    let mut lattice = SimulationLattice::new(LATTICE_SIZE);

    // Store all frames in memory
    let mut frames = Vec::new();

    for step in 0..TIME_STEPS {
        lattice.evolve(step);
        let frame = lattice.create_2d_flatmap();
        frames.push(frame);
    }

    // Consolidate all frames into one PNG
    // Here we stack them vertically in one column
    let single_frame_width = frames[0].width();
    let single_frame_height = frames[0].height();
    let total_height = single_frame_height * (frames.len() as u32);
    let total_width = single_frame_width;

    let mut big_image = RgbImage::new(total_width, total_height);
    for (i, frame) in frames.iter().enumerate() {
        let y_offset = i as u32 * single_frame_height;
        for x in 0..single_frame_width {
            for y in 0..single_frame_height {
                let pixel = frame.get_pixel(x, y);
                big_image.put_pixel(x, y_offset + y, *pixel);
            }
        }
    }

    big_image.save("all_frames_consolidated.png").unwrap();

    // Export 3D data
    lattice.export_3d_data("3d_data.csv");

    println!("Simulation complete. All frames consolidated into all_frames_consolidated.png");
    println!("3D data exported as 3d_data.csv.");
}
