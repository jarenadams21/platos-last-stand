use std::f64::consts::PI;

/// Constants
const HBAR: f64 = 6.582119569e-16; // Reduced Planck constant in eV·s
const C: f64 = 299_792_458.0;       // Speed of light in m/s
const ELECTRON_MASS: f64 = 0.5109989461e6; // Electron mass in eV/c^2

/// Enum to represent particle types
#[derive(Debug, Clone)]
enum ParticleType {
    Fermion,
    Boson,
}

/// Struct to represent a quantum state
#[derive(Debug, Clone)]
struct QuantumState {
    particle_type: ParticleType,
    mass: f64,    // in eV/c^2
    energy: f64,  // in eV
    momentum: f64, // in eV·s/m
    position: f64, // in meters
    spin: f64,
}

impl QuantumState {
    /// Initialize a new quantum state
    fn new(
        particle_type: ParticleType,
        mass: f64,
        energy: f64,
        momentum: f64,
        position: f64,
        spin: f64,
    ) -> Self {
        QuantumState {
            particle_type,
            mass,
            energy,
            momentum,
            position,
            spin,
        }
    }

    /// Calculate the Compton wavelength
    fn compton_wavelength(&self) -> f64 {
        HBAR / (self.mass * C)
    }
}

/// Operators for creation and annihilation
#[derive(Debug)]
struct Operators;

impl Operators {
    /// Annihilation operator
    fn annihilate(state: &QuantumState) -> Option<QuantumState> {
        // In a full quantum field theory, this would reduce the particle number
        // For simulation, we'll represent this as returning None
        None
    }

    /// Creation operator
    fn create(particle_type: ParticleType, mass: f64, spin: f64) -> QuantumState {
        // Creates a new particle with given properties
        let energy = mass * C.powi(2);
        let momentum = 0.0; // At rest
        let position = 0.0; // Arbitrary
        QuantumState::new(particle_type, mass, energy, momentum, position, spin)
    }
}

/// Hamiltonian operator representing the total energy
struct Hamiltonian;

impl Hamiltonian {
    /// Applies the Hamiltonian to evolve the state over time
    fn evolve(state: &QuantumState, time: f64) -> QuantumState {
        // Time evolution according to Schrödinger's equation (simplified)
        // For stationary states, the state remains unchanged
        state.clone()
    }
}

/// System representing the entire quantum system
struct System {
    states: Vec<QuantumState>,
}

impl System {
    /// Initialize a new system
    fn new() -> Self {
        System { states: Vec::new() }
    }

    /// Add a quantum state to the system
    fn add_state(&mut self, state: QuantumState) {
        self.states.push(state);
    }

    /// Simulate annihilation of an electron and positron
    fn simulate_annihilation(&mut self) {
        // Find electron and positron indices
        let mut electron_index = None;
        let mut positron_index = None;

        for (i, state) in self.states.iter().enumerate() {
            match state.particle_type {
                ParticleType::Fermion => {
                    if state.mass == ELECTRON_MASS && state.spin == 0.5 {
                        if state.energy > 0.0 {
                            electron_index = Some(i);
                        } else {
                            positron_index = Some(i);
                        }
                    }
                }
                _ => {}
            }
        }

        // Proceed if both are found
        if let (Some(e_idx), Some(p_idx)) = (electron_index, positron_index) {
            // Remove electron and positron from states
            // Remove higher index first to avoid shifting
            if e_idx > p_idx {
                self.states.remove(e_idx);
                self.states.remove(p_idx);
            } else {
                self.states.remove(p_idx);
                self.states.remove(e_idx);
            }

            // Apply annihilation operator
            let annihilation_result = Operators::annihilate(&QuantumState::new(
                ParticleType::Fermion,
                ELECTRON_MASS,
                0.0,
                0.0,
                0.0,
                0.5,
            ));

            // Simulate creation of photons (massless bosons)
            let photon_energy = 2.0 * ELECTRON_MASS * C.powi(2); // Total energy converted
            let photon = Operators::create(ParticleType::Boson, 0.0, 1.0);
            let mut photon1 = photon.clone();
            let mut photon2 = photon.clone();
            photon1.energy = photon_energy / 2.0;
            photon2.energy = photon_energy / 2.0;

            // Add photons to the system
            self.states.push(photon1);
            self.states.push(photon2);
        } else {
            println!("Electron and positron not found for annihilation.");
        }
    }

    /// Log the current state of the system
    fn log_system_state(&self) {
        println!("System State:");
        for state in &self.states {
            println!("{:?}", state);
        }
    }
}

fn main() {
    // Initialize the system
    let mut system = System::new();

    // Create an electron and a positron
    let electron = QuantumState::new(
        ParticleType::Fermion,
        ELECTRON_MASS,
        ELECTRON_MASS * C.powi(2),
        0.0,
        0.0,
        0.5,
    );

    let positron = QuantumState::new(
        ParticleType::Fermion,
        ELECTRON_MASS,
        -ELECTRON_MASS * C.powi(2),
        0.0,
        0.0,
        0.5,
    );

    let electron2 = QuantumState::new(
        ParticleType::Fermion,
        ELECTRON_MASS,
        ELECTRON_MASS * C.powi(2),
        0.0,
        0.0,
        0.5,
    );

    let electron3 = QuantumState::new(
        ParticleType::Fermion,
        ELECTRON_MASS,
        ELECTRON_MASS * C.powi(2),
        0.0,
        0.0,
        0.5,
    );

    let electron4 = QuantumState::new(
        ParticleType::Fermion,
        ELECTRON_MASS,
        ELECTRON_MASS * C.powi(2),
        0.0,
        0.0,
        0.5,
    );

    // Add them to the system
    system.add_state(electron);
    system.add_state(positron);
    system.add_state(electron2);
    system.add_state(electron3);
    system.add_state(electron4);

    // Log initial state
    println!("Initial System State:");
    system.log_system_state();

    // Simulate annihilation
    system.simulate_annihilation();

    // Log final state
    println!("\nFinal System State after Annihilation:");
    system.log_system_state();
}
