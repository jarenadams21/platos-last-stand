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
    pub total_energy: f64,       // in MeV
    pub total_mass: f64,         // in MeV/c^2
    pub dark_energy_percentage: f64,
    pub dark_matter_percentage: f64,
    pub atoms_percentage: f64,
}

impl System {
    /// Creates a new system.
    pub fn new() -> Self {
        System {
            particles: Vec::new(),
            total_energy: 0.0,
            total_mass: 0.0,
            dark_energy_percentage: 0.0,
            dark_matter_percentage: 0.0,
            atoms_percentage: 0.0,
        }
    }

    /// Adds a particle to the system.
    pub fn add_particle(&mut self, particle: Particle) {
        self.total_mass += particle.mass;
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
            // Get the masses before removal
            let electron_mass = self.particles[e_index].mass;
            let positron_mass = self.particles[p_index].mass;

            // Remove electron and positron
            // Remove the higher index first to avoid index shifting
            if e_index > p_index {
                self.particles.remove(e_index);
                self.particles.remove(p_index);
            } else {
                self.particles.remove(p_index);
                self.particles.remove(e_index);
            }

            // Total rest mass energy released
            let total_mass_energy = electron_mass + positron_mass; // in MeV

            // Create massless long-range binding force particle (Graviton)
            let graviton = Particle::new(
                "Graviton",
                ParticleType::Boson,
                0.0,     // Massless
                0.0,     // Charge
                2.0,     // Spin 2
                1,       // Point-split density of 1 (massless and chargeless)
            );

            // Create massless energy momentum boson (Photon)
            let photon = Particle::new(
                "Photon",
                ParticleType::Boson,
                0.0,     // Massless
                0.0,
                1.0,
                1,       // Point-split density of 1 (massless and chargeless)
            );

            self.particles.push(graviton);
            self.particles.push(photon);

            // Update total mass and energy
            self.total_mass -= total_mass_energy;
            self.total_energy += total_mass_energy;

            // Simulate emergent particles from graviton and photon interaction
            self.simulate_emergent_particles(total_mass_energy);

            // Update percentages based on energy content
            self.update_percentages();

            // Log observables
            // self.log_observables();
        } else {
            println!("Electron and positron not found in the system.");
        }
    }

    /// Simulates the interaction of graviton and photon producing emergent particles.
    pub fn simulate_emergent_particles(&mut self, available_energy: f64) {
        // For simplicity, assume that the available energy can be used to create new particles
        // Here, we'll create an electron and a neutrino

        let electron_mass = 0.511; // MeV/c^2
        let neutrino_mass = 0.0;   // Approximate neutrino mass as negligible

        // Check if enough energy is available to create an electron
        if available_energy >= electron_mass {
            // Create an electron
            let electron = Particle::new(
                "Electron",
                ParticleType::Fermion,
                electron_mass,
                -1.0,
                0.5,
                3,       // Point-split density of 3
            );
            self.particles.push(electron);

            // Update total mass and energy
            self.total_mass += electron_mass;
            self.total_energy -= electron_mass; // Energy used to create mass

            // Create a neutrino with negligible mass
            let neutrino = Particle::new(
                "Neutrino",
                ParticleType::Fermion,
                neutrino_mass,
                0.0,
                0.5,
                1,       // Point-split density of 1
            );
            self.particles.push(neutrino);

            // Remaining energy is carried away as kinetic energy of particles
            // For simplicity, we won't track kinetic energy here

        } else {
            println!("Not enough energy to create emergent particles.");
        }
    }

    /// Updates the system's energy and matter percentages.
    pub fn update_percentages(&mut self) {
        // Calculate total mass-energy content
        let total_mass_energy = self.particles.iter().fold(0.0, |sum, p| sum + p.mass);

        // For simplicity, assume:
        // - Dark energy is associated with the energy released (massless particles and energy in the system).
        // - Dark matter remains constant (since no dark matter particles are involved).
        // - Atomic matter is associated with the massive particles.

        let dark_energy = self.total_energy;
        let atoms_mass_energy = total_mass_energy;
        let total_content = dark_energy + atoms_mass_energy;

        self.dark_energy_percentage = (dark_energy / total_content) * 100.0;
        self.atoms_percentage = (atoms_mass_energy / total_content) * 100.0;

        // Assuming dark matter percentage remains the same as before the annihilation
        // For simplicity, we can set it to balance to 100%
        self.dark_matter_percentage = 100.0 - self.dark_energy_percentage - self.atoms_percentage;
    }

    /// Logs the current observables and particle states.
    pub fn log_observables(&self) {
        println!("System Observables:");
        println!("Total Energy Released: {:.3} MeV", self.total_energy);
        println!("Total Mass Remaining: {:.3} MeV/c^2", self.total_mass);
        println!("Dark Energy: {:.3}%", self.dark_energy_percentage);
        println!("Dark Matter: {:.3}%", self.dark_matter_percentage);
        println!("Atoms: {:.3}%", self.atoms_percentage);
        println!("Particles in the system:");
        for particle in &self.particles {
            println!("{:?}", particle);
        }
    }
}
