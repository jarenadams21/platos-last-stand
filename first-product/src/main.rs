use ndarray::prelude::*;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use lazy_static::lazy_static;

// Fundamental constants matching 1:1 with the Standard Model and relativistic units
// Units: c = ħ = k_B = 1, M_P = (8πG)^(-1/2)
const MP: f64 = 2.435e18; // Reduced Planck mass in GeV
const M: f64 = 1.0e9;     // Example scale, adjustable to tested scenario
const HI: f64 = MP;       // Primeval scale at or below reduced Planck mass
const H0: f64 = 1.0e-33;  // Present-day Hubble scale (~ in GeV)
const NU: f64 = 0.001;    // Example small deviation parameter from rigid vacuum
const OMEGA_M0: f64 = 0.3;
const OMEGA_R0: f64 = 1e-5;

lazy_static! {
    static ref HF: f64 = H0*((OMEGA_M0/(1.0 - NU)).sqrt()); // Final de Sitter scale derived
}

// Equation of state indices for radiation and matter
// Radiation: ω = 1/3, Matter: ω = 0
// The code runs through radiation-to-matter and matter-to-final de Sitter phases.

// Derived parameters from the theoretical appendices:
fn alpha_param() -> f64 {
    // According to the final expression from the appendices:
    // α = ((1−ν)^3(1−2ν))^{1/4} * (HI/(2^{1/4} * MP * M))
    let one_minus_nu = 1.0 - NU;
    let one_minus_2nu = 1.0 - 2.0*NU;
    let val = ((one_minus_nu.powf(3.0)) * one_minus_2nu).powf(0.25);
    val * (HI/(M*(2.0_f64.powf(0.25)*MP)))
}

fn sigma_param() -> f64 {
    // For matter era (ω=0):
    // σ = (1−ν)*HF/(M^2)
    (1.0 - NU)*(*HF)/(M*M)
}

// Potential and kinetic energy densities in the radiation era:
fn potential_radiation(phi: f64, alpha: f64) -> f64 {
    // V(φ)=VI [1+(α φ)^4]^-2 [1+ν(α φ)^4]
    let vi = 3.0*MP.powi(2)*HI.powi(2);
    let x = (alpha*phi).powi(4);
    vi * (1.0 + NU*x)/((1.0 + x).powi(2))
}

fn kinetic_radiation(phi: f64, alpha: f64) -> f64 {
    // ρ_k(φ) = (1−ν)*VI*(αφ)^4 / [ (1+(αφ)^4)^2 ]
    let vi = 3.0*MP.powi(2)*HI.powi(2);
    let x = (alpha*phi).powi(4);
    (1.0 - NU)*vi*x/((1.0 + x).powi(2))
}

// Potential and kinetic energy densities in the matter era:
fn potential_matter(phi: f64, sigma: f64) -> f64 {
    // V(φ)=VF [1+ ν(σ φ)^{-3}]
    // VF = 3 MP^2 HF^2
    let vf = 3.0*MP.powi(2)*(*HF).powi(2);
    let y = (sigma*phi).powf(-3.0);
    vf*(1.0 + NU*y)
}

fn kinetic_matter(phi: f64, sigma: f64) -> f64 {
    // ρ_k(φ) = (1−ν)*VF (σ φ)^{-3}/[1+ν(σ φ)^{-3}]
    let vf = 3.0*MP.powi(2)*(*HF).powi(2);
    let y = (sigma*phi).powf(-3.0);
    (1.0 - NU)*vf*y/(1.0 + NU*y)
}

// Time evolution parameters:
const TIME_STEPS: usize = 1000;
const DT: f64 = 5.59 * 1e-44; // Very small timestep to resolve early universe dynamics

struct FieldPoint {
    phi: f64,
    phidot: f64,
}

// Evolve during radiation era:
fn evolve_radiation(phi: f64, phidot: f64, alpha: f64) -> (f64, f64) {
    // H(φ)=HI/[1+(αφ)^4]^{1/2}
    let x = (alpha*phi).powi(4);
    let h_val = HI/(1.0 + x).sqrt();

    // From eq.(56): φ̇ = (1−ν)*φ*H(φ)
    let phidot_new = (1.0 - NU)*phi*h_val;

    // To find φ̈, we differentiate φ̇ w.r.t φ:
    // φ̇(φ) = (1−ν)*H(φ)*φ
    // dφ̇/dφ = (1−ν)(H(φ) + φ dH/dφ)
    // dH/dφ = HI*(-2)*α^4 φ^3/(1+x)^{3/2}
    let alpha4 = alpha.powi(4);
    let dH_dphi = HI*(-2.0)*alpha4*phi.powi(3)/((1.0+x).powf(1.5));
    let dphidot_dphi = (1.0 - NU)*(h_val + phi*dH_dphi);

    // φ̈ = dφ̇/dt = dφ̇/dφ * φ̇
    let phiddot = dphidot_dphi * phidot_new;

    let phi_next = phi + phidot*DT;
    let phidot_next = phidot + phiddot*DT;

    (phi_next, phidot_next)
}

// Evolve during matter era:
fn evolve_matter(phi: f64, phidot: f64, sigma: f64) -> (f64, f64) {
    // H(φ)=HF(1+(σφ)^{-3})^{1/2}
    let y = (sigma*phi).powf(-3.0);
    let h_val = (*HF)*(1.0+y).sqrt();

    let phidot_new = (1.0 - NU)*phi*h_val;

    // dH/dφ for matter era:
    // y=(σ φ)^{-3}, dy/dφ = -3σ^{-3}φ^{-4}
    let dy_dphi = -3.0*(sigma.powf(-3.0))*phi.powf(-4.0);
    let dH_dphi = (*HF)*0.5*(1.0+y).powf(-0.5)*dy_dphi;

    let dphidot_dphi = (1.0 - NU)*(h_val + phi*dH_dphi);
    let phiddot = dphidot_dphi * phidot_new;

    let phi_next = phi + phidot*DT;
    let phidot_next = phidot + phiddot*DT;

    (phi_next, phidot_next)
}

fn main() {
    let alpha = alpha_param();
    let sigma = sigma_param();

    // Initial conditions:
    let mut phi = 77.3147 * (1.0/137.0);
    let x_init = (alpha*1e-30).powi(4);
    let h_initial = HI/(1.0+x_init).sqrt();
    let mut phidot = (1.0 - NU)*1e-30*h_initial;

    let phi_end_radiation = 1.0/alpha;

    let mut file = File::create("qgp_data.csv").unwrap();
    writeln!(file, "step,phi,phidot,H,Potential,Kinetic,Era").unwrap();

    // Radiation era
    for step in 0..TIME_STEPS {
        let pot = potential_radiation(phi, alpha);
        let kin = kinetic_radiation(phi, alpha);
        let x = (alpha*phi).powi(4);
        let h_val = HI/(1.0+x).sqrt();
        writeln!(file, "{},{},{},{},{},{},{}",
                 step, phi, phidot, h_val, pot, kin, "radiation").unwrap();

        if phi > phi_end_radiation {
            break;
        }
        let (phi_new, phidot_new) = evolve_radiation(phi, phidot, alpha);
        phi = phi_new;
        phidot = phidot_new;
    }

    // Matter era
    for step in TIME_STEPS..(2*TIME_STEPS) {
        let pot = potential_matter(phi, sigma);
        let kin = kinetic_matter(phi, sigma);
        let y = (sigma*phi).powf(-3.0);
        let h_val = (*HF)*(1.0+y).sqrt();

        writeln!(file, "{},{},{},{},{},{},{}",
                 step, phi, phidot, h_val, pot, kin, "matter").unwrap();

        let (phi_new, phidot_new) = evolve_matter(phi, phidot, sigma);
        phi = phi_new;
        phidot = phidot_new;

        if (h_val - *HF).abs() < 1e-35 {
            break;
        }
    }

    println!("Simulation complete. Data in quantum_gravity_plasma_data.csv");
}
