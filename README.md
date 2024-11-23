# Quantum Annihilation Simulation Using Heisenberg-Lee Model in Rust

## Introduction

This project simulates a quantum annihilation experiment using **Heisenberg's non-linear spinor theory** within the framework of the **Lee model**. It models the collision and annihilation of an electron and a positron, resulting in the creation of massless photons. The simulation tracks observables such as dark energy, dark matter, and atomic matter, conforming to the latest Planck statistics.

Implemented in Rust, this codebase provides a flexible and extensible foundation for further development and experimentation in quantum field theory simulations.

---

## Code Correspondence

The theoretical framework is implemented in Rust as follows:

### Particle Representation

- **`Particle` Struct**:

  - Represents particles with attributes like name, type (fermion or boson), mass, charge, spin, point-split density, and length dimension.
  - **Length Dimension Calculation**: Based on point-split density, adhering to the classification:

    ```rust
    let length_dimension = match point_split_density {
        0 | 1 => -3.0 / 2.0, // Psi_3 objects
        2 => -1.0,           // Psi_2 objects
        3 => -1.5,           // Psi_3 objects
        4 => -2.0,           // Psi_4 objects
        _ => -2.0,           // Unstable particles (Psi_4)
    };
    ```

### System Simulation

- **`System` Struct**:

  - Manages the collection of particles and observables.
  - **`simulate_annihilation` Method**:

    - Simulates the annihilation of an electron and positron into two photons.
    - Models the interaction Hamiltonian \( H_{\text{int}} \).

    ```rust
    pub fn simulate_annihilation(&mut self) {
        // ... find and remove electron and positron
        // ... create two photons (gamma particles)
        // ... update observables
    }
    ```

### Observables Tracking

- **Dark Energy, Dark Matter, Atoms Percentages**:

  - Updated to match Planck statistics after the annihilation event:

    ```rust
    self.dark_energy_percentage = 68.3;
    self.dark_matter_percentage = 26.8;
    self.atoms_percentage = 4.9;
    ```

- **Logging**:

  - The `log_observables` method outputs the current state of the system, including particle details and observables.

### Stability Checks

- **`is_stable` Method**:

  - Determines particle stability based on point-split density:

    ```rust
    pub fn is_stable(&self) -> bool {
        self.point_split_density <= 4
    }
    ```

### Major Laws Implementation

- **Conservation Laws**:

  - Ensured by the logic in particle interactions and annihilations.

- **Pauli Exclusion Principle**:

  - Reflected in the handling of fermionic particles, avoiding duplicate states.

- **Symmetry Principles**:

  - The code structure allows for symmetric operations and can be extended to include PT symmetry considerations.

---

## Further Development Guidelines

To facilitate easy extension and enhancement of the simulation, consider the following:

### 1. Extending Particle Types

- **Add New Particles**:

  - Introduce additional particles by creating new instances of the `Particle` struct.
  - Ensure point-split density and length dimensions are correctly assigned.

### 2. Implementing Additional Interactions

- **Define New Interaction Methods**:

  - Create methods within the `System` struct to simulate other particle interactions (e.g., decay processes, scattering).

- **Respect Conservation Laws**:

  - Ensure all new interactions conserve energy, momentum, and charge.

### 3. Enhancing Observables Tracking

- **Detailed Metrics**:

  - Implement tracking of additional observables like entropy, temperature, or particle lifetimes.

- **Visualization**:

  - Integrate visualization tools to represent the system's evolution graphically.

### 4. Incorporating Advanced Theoretical Concepts

- **Quantum Field Theory Elements**:

  - Expand the simulation to include concepts like spontaneous symmetry breaking, gauge fields, and renormalization.

- **Non-Hermitian Hamiltonians and PT Symmetry**:

  - Explore regimes where the Hamiltonian is non-Hermitian and analyze the implications of PT symmetry in the simulation.

### 5. Code Optimization and Testing

- **Performance Optimization**:

  - Profile the code to identify bottlenecks and optimize accordingly.

- **Testing Framework**:

  - Develop unit tests and integration tests to ensure the correctness of simulations.

- **Documentation**:

  - Maintain thorough documentation for all code components to aid in future development.

---

## References

- **Heisenberg, W.** *Non-Linear Spinor Theory*. (Original papers on non-linear spinor theory).
- **Lee, T.D.** *Particle Physics and Introduction to Field Theory*. (For insights into the Lee model).
- **Planck Collaboration**. *Planck 2018 Results*. (For the latest cosmological observations and statistics).

---

**Note**: This simulation is a simplified representation and serves as a foundation for educational and exploratory purposes in quantum physics simulations. :D