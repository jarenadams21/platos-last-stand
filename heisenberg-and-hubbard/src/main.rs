// main.rs

mod lib;

use lib::{Particle, ParticleType, System};

fn main() {
    let mut system = System::new();

    // Create initial particles: electron, positron, and neutrino
    let electron = Particle::new("Electron", ParticleType::Fermion, 0.511, -1.0, 0.5, 3);
    let positron = Particle::new("Positron", ParticleType::Fermion, 0.511, 1.0, 0.5, 3);
    let neutrino = Particle::new("Neutrino", ParticleType::Fermion, 0.0, 0.0, 0.5, 1);

    // Add particles to the system
    system.add_particle(electron);
    system.add_particle(positron);
    system.add_particle(neutrino);

    println!("Before annihilation:");
    system.log_observables();

    // Simulate the annihilation
    system.simulate_annihilation();

    println!("\nAfter annihilation:");
    system.log_observables();
}
