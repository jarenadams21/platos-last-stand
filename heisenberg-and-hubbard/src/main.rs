// main.rs

mod lib;

use lib::{Particle, ParticleType, System};

fn main() {
    let mut system = System::new();

    // Create initial particles: electron and positron
    let electron = Particle::new("Electron", ParticleType::Fermion, 0.511, -1.0, 0.5, 3);
    let positron = Particle::new("Positron", ParticleType::Fermion, 0.511, 1.0, 0.5, 3);

    // Add particles to the system
    system.add_particle(electron);
    system.add_particle(positron);

    println!("Before annihilation:\n\n");
    system.log_observables();

    // Simulate the annihilation
    system.simulate_annihilation();

    println!("\nAfter annihilation:\n\n");
    system.log_observables();
}

/*

mod lib;

use lib::{Particle, ParticleType, System};

fn main() {
    let mut system = System::new();

    // Create initial particles: electron, positron, neutrino, and a proton
    let electron = Particle::new("Electron", ParticleType::Fermion, 0.511, -1.0, 0.5, 3);
    let positron = Particle::new("Positron", ParticleType::Fermion, 0.511, 1.0, 0.5, 3);
    let neutrino = Particle::new("Neutrino", ParticleType::Fermion, 0.0, 0.0, 0.5, 1);
    let proton = Particle::new("Proton", ParticleType::Fermion, 938.272, 1.0, 0.5, 3);

    // Add particles to the system
    system.add_particle(electron);
    system.add_particle(positron);
    system.add_particle(neutrino);
    system.add_particle(proton);

    println!("Before annihilation:");
    system.log_observables();

    // Simulate the annihilation
    system.simulate_annihilation();

    println!("\nAfter annihilation:");
    system.log_observables();
}
*/