// lib.rs

use std::f64::consts::PI;

/// Enum to represent the type of particle: Fermion or Boson.
#[derive(Debug, Clone)]
pub enum ParticleType {
    Fermion,
    Boson,
}

/// Struct to represent a particle with its properties.
#[derive(Debug, Clone)]
pub struct Particle {
    pub name: String,
    pub particle_type: ParticleType,
    pub mass: f64,               // in MeV/c^2
    pub charge: f64,             // in elementary charge units
    pub spin: f64,               // spin quantum number
    pub point_split_density: u32,
    pub length_dimension: f64,   // Psi_n where n corresponds to length dimension
}

impl Particle {
    /// Creates a new particle with calculated length dimension.
    pub fn new(
        name: &str,
        particle_type: ParticleType,
        mass: f64,
        charge: f64,
        spin: f64,
        point_split_density: u32,
    ) -> Self {
        let length_dimension = match point_split_density {
            0 | 1 => -3.0 / 2.0,     // Psi_3 objects
            2 => -1.0,               // Psi_2 objects
            3 => -1.5,               // Psi_3 objects
            4 => -2.0,               // Psi_4 objects
            _ => -2.0,               // Unstable particles (Psi_4)
        };

        Particle {
            name: name.to_string(),
            particle_type,
            mass,
            charge,
            spin,
            point_split_density,
            length_dimension,
        }
    }

    /// Determines if a particle is stable based on point-split density.
    pub fn is_stable(&self) -> bool {
        self.point_split_density <= 4
    }
}

/// Struct to represent the system containing particles and observables.
#[derive(Debug)]
pub struct System {
    pub particles: Vec<Particle>,
    pub dark_energy_percentage: f64,
    pub dark_matter_percentage: f64,
    pub atoms_percentage: f64,
}

impl System {
    /// Creates a new system.
    pub fn new() -> Self {
        System {
            particles: Vec::new(),
            dark_energy_percentage: 0.0,
            dark_matter_percentage: 0.0,
            atoms_percentage: 0.0,
        }
    }

    /// Adds a particle to the system.
    pub fn add_particle(&mut self, particle: Particle) {
        self.particles.push(particle);
    }

    /// Simulates the annihilation of an electron and positron.
    pub fn simulate_annihilation(&mut self) {
        // Find electron and positron indices
        let mut electron_index = None;
        let mut positron_index = None;

        for (i, particle) in self.particles.iter().enumerate() {
            if particle.name == "Electron" && particle.charge == -1.0 {
                electron_index = Some(i);
            } else if particle.name == "Positron" && particle.charge == 1.0 {
                positron_index = Some(i);
            }
        }

        // Proceed if both particles are found
        if let (Some(e_index), Some(p_index)) = (electron_index, positron_index) {
            // Remove electron and positron
            self.particles.remove(e_index);
            // Adjust index if positron was after electron
            if p_index > e_index {
                self.particles.remove(p_index - 1);
            } else {
                self.particles.remove(p_index);
            }

            // Create two photons (gamma particles)
            let gamma1 = Particle::new(
                "Photon",
                ParticleType::Boson,
                0.0,
                0.0,
                1.0,
                1, // Point-split density of 1 (massless and chargeless)
            );

            let gamma2 = gamma1.clone();

            self.particles.push(gamma1);
            self.particles.push(gamma2);

            // Update percentages based on Planck statistics
            self.update_percentages();

            // Log observables
            self.log_observables();
        } else {
            println!("Electron and positron not found in the system.");
        }
    }

    /// Updates the system's energy and matter percentages.
    pub fn update_percentages(&mut self) {
        // Updating percentages to match Planck statistics
        self.dark_energy_percentage = 68.3;
        self.dark_matter_percentage = 26.8;
        self.atoms_percentage = 4.9;
    }

    /// Logs the current observables and particle states.
    pub fn log_observables(&self) {
        println!("System Observables:");
        println!("Dark Energy: {}%", self.dark_energy_percentage);
        println!("Dark Matter: {}%", self.dark_matter_percentage);
        println!("Atoms: {}%", self.atoms_percentage);
        println!("Particles in the system:");
        for particle in &self.particles {
            println!("{:?}", particle);
        }
    }
}
