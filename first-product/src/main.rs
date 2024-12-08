use ndarray::prelude::*;
use rand::SeedableRng;
use rand::{rngs::StdRng, Rng};
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use rug::Float; // For high-precision arithmetic

// Fundamental constants matching 1:1 with the Standard Model and relativistic units
// Units: c = ħ = k_B = 1, M_P = (8πG)^(-1/2)
const MP: f64 = 2.435e18; // Reduced Planck mass in GeV
const M: f64 = 1.0e9;     // Example scale, adjustable to tested scenario
const HI: f64 = MP;       // Primeval scale at or below reduced Planck mass
const H0: f64 = 1.0e-33;  // Present-day Hubble scale (~ in GeV)
const NU: f64 = 0.001;    // Example small deviation parameter from rigid vacuum
const OMEGA_M0: f64 = 0.3;
const OMEGA_R0: f64 = 1e-5;
const HF: f64 = H0*((OMEGA_M0/(1.0 - NU)).sqrt()); // Final de Sitter scale derived

// Equation of state indices for radiation and matter
// Radiation: ω = 1/3, Matter: ω = 0
// The code runs through radiation-to-matter and matter-to-final de Sitter phases.

// Non-canonical scalar field parameters for the radiation era (β=2):
// Potential: V(φ)=VI*[1+(αφ)^4]^-2[1+ν(αφ)^4]
// where VI = 3 * MP^2 * HI^2
// α determined as derived in Appendix C
fn alpha_param() -> f64 {
    // As derived: α = [2/( (1-ν)^3(1-2ν) )^(1/4)] * (M/MP)*(HI^-1)
    // The derived expression from the paper:
    // α = ( (1−ν)^3 (1−2ν) )^(1/4) * sqrt( (2 * M * MP) / HI )
    // Adjusting from the provided detailed appendices:
    // Actually final formula from code snippet in the text:
    // α= (1/(2^(1/4))) * ((1−ν)^(3/4) (1−2ν)^(1/4)) * (MP/M)*(1/HI)
    // There was a slight misalignment, but we use the exact final from appendix C:
    let one_minus_nu = 1.0 - NU;
    let one_minus_2nu = 1.0 - 2.0*NU;
    // The text states:
    // φe = α^-1 = [2 (1−ν)^3(1−2ν)]^(1/4) * (M/ (sqrt(1−ν)*MP*...)) etc.
    // We follow the final given formula strictly (C7):
    // α = [ ( (1−ν)^3(1−2ν) )^(1/4) * sqrt(MP) * (2^(−1/4)) ] / (M * (HI/M) )
    // Reviewing the paper steps carefully, final consistent choice:
    let val = ((one_minus_nu).powf(3.0)*(one_minus_2nu)).powf(0.25);
    // From the paper: α = [ (1−ν)^{3}(1−2ν) ]^{1/4} * (HI/M) / (2^(1/4)*MP)
    // Actually it gave a factor in a slightly different form. We must ensure consistency.
    // We want sub-Planckian: Let's trust the final given expression (C7):
    // α= ((1−ν)^{3}(1−2ν))^{1/4} * (HI/(2^(1/4)*MP*M))
    let alpha = val * (HI/(M*((2.0_f64).powf(0.25)*MP)));
    alpha
}

// For the matter era (ω=0), we use σ parameter from (C12):
fn sigma_param() -> f64 {
    // σ = (1−ν)*HF/(M^2), from the text in the limit ω->0
    (1.0 - NU)*HF/(M*M)
}

// Compute potential and kinetic energy densities for given φ in radiation phase:
fn potential_radiation(phi: f64, alpha: f64) -> f64 {
    // V(φ)=VI [1+(α φ)^4]^-2 [1+ν(α φ)^4]
    let vi = 3.0*MP.powi(2)*HI.powi(2);
    let x = (alpha*phi).powi(4);
    vi * (1.0 + NU*x)/((1.0 + x).powi(2))
}
fn kinetic_radiation(phi: f64, alpha: f64) -> f64 {
    // ρ_k(φ) = (1−ν) VI (αφ)^4 / [ (1+(αφ)^4)^2 ]
    let vi = 3.0*MP.powi(2)*HI.powi(2);
    let x = (alpha*phi).powi(4);
    (1.0 - NU)*vi*(x)/((1.0 + x).powi(2))
}

// In the matter era:
fn potential_matter(phi: f64, sigma: f64) -> f64 {
    // V(φ)=VF [1+ ν(σ φ)^{-3}]
    // VF = 3 MP^2 HF^2
    let vf = 3.0*MP.powi(2)*HF.powi(2);
    let y = (sigma*phi).powf(-3.0);
    vf*(1.0 + NU*y)
}
fn kinetic_matter(phi: f64, sigma: f64) -> f64 {
    // ρ_k(φ) = (1−ν)*VF (σ φ)^{-3}/[1+ν(σ φ)^{-3}]
    let vf = 3.0*MP.powi(2)*HF.powi(2);
    let y = (sigma*phi).powf(-3.0);
    (1.0 - NU)*vf*y/(1.0 + NU*y)
}

// Time evolution parameters
const TIME_STEPS: usize = 1000;
const DT: f64 = 1e-35; // Example time step, very small to resolve Planck-scale dynamics

// We define a structure for our fields:
#[derive(Clone, Debug)]
struct FieldPoint {
    phi: Float,
    phidot: Float,
    // Additional fields as needed for full quantum gravity plasma states:
    // We consider radiation and matter densities as derived from φ.
}

// We define a function to update φ based on the equations of motion:
// For radiation era (early times), we have φ evolution from equations (C3),(C4) etc.
// For matter era (later times), we switch the ω accordingly.
// The code will run from HI scale down to HF scale, using the derived equations.

fn evolve_radiation(phi: &Float, phidot: &Float, alpha: f64) -> (Float, Float) {
    // Equation of motion from section V and appendices for radiation era.
    // φ̇ and φ̈ derived from conditions in the text. Here we numerically approximate:
    let phi_val = phi.clone();
    let phidot_val = phidot.clone();

    // For simplicity, we integrate a known solution path. The text provides a closed form for H(φ):
    // H(φ)= (HI)/[1+(α φ)^4]^{1/2}
    // From eq. (56)
    // To find φ̈: differentiate relations numerically or use a derived formula from references.

    // Since full closed-form φ(t) is complex, we will approximate with finite differences or a chosen ODE method.
    // We'll pick a midpoint method to ensure stability:

    let alpha_f = Float::with_val(128, alpha);
    let hi_f = Float::with_val(128, HI);
    let nu_f = Float::with_val(128, NU);
    let one = Float::with_val(128, 1.0);

    let x = (&alpha_f * &phi_val).pow(4);
    // H(φ) in radiation era:
    let h_val = &hi_f/(one.clone() + &x).sqrt();

    // For the equation of motion, we have:
    // φ̇ = (from eq. (56), φ̇ = (1−ν)φ H(φ))
    // Actually from eq.(56): H(φ) = φ̇/[(1−ν)φ], so φ̇ = H(φ)*(1−ν)*φ
    let one_minus_nu = Float::with_val(128, 1.0 - NU);
    let phidot_new = &h_val * &one_minus_nu * &phi_val;

    // φ̈ requires differentiating phidot_new wrt φ and chain rule with dφ/dt.
    // Numerically approximate φ̈ = dφ̇/dt = dφ̇/dφ * dφ/dt.
    // dφ̇/dφ = derivative of [H(φ)*(1−ν)*φ]
    // H(φ) = HI/(1+x)^{1/2}, x=(α φ)^4
    // dH/dφ = HI*(-(1/2)*(1+x)^{-3/2}*(4(α^4)φ^3))
    let alpha4 = alpha_f.pow(4);
    let dH_dphi = &hi_f * Float::with_val(128, -2.0)*&alpha4*phi_val.pow(3)/(Float::with_val(128,2.0)*(one.clone()+&x).pow(3.0/2.0));
    let dH_dphi_corrected = &dH_dphi * Float::with_val(128, -2.0); 
    // Correction: Check chain rule carefully:
    // dH/dφ = d/dφ [HI(1+x)^{-1/2}] with x=(αφ)^4
    // dx/dφ=4α^4 φ^3
    // dH/dφ = HI * (-1/2)*(1+x)^(-3/2)*4α^4 φ^3
    // dH/dφ = HI*(-2)*α^4 φ^3 /(1+x)^{3/2}
    let dH_dphi_final = Float::with_val(128, -2.0)*&hi_f*&alpha4*phi_val.pow(3)/(one.clone()+&x).pow(3.0/2.0);

    let dphidot_dphi = &one_minus_nu*(&h_val + &phi_val*&dH_dphi_final);
    let phiddot = &dphidot_dphi * &phidot_new;

    let phi_next = phi_val + &phidot_val*Float::with_val(128, DT);
    let phidot_next = phidot_val + phiddot*Float::with_val(128, DT);

    (phi_next, phidot_next)
}

fn evolve_matter(phi: &Float, phidot: &Float, sigma: f64) -> (Float, Float) {
    // Similar approach for matter era (ω=0).
    // H(φ)= HF[1+D(1+z)^{3(1−ν)}]^{1/2} but we directly use eq.62 form:
    // From eq.(62): H(φ) = φ̇/[(1−ν)φ], so again φ̇ = (1−ν)*φ*H(φ).
    // For matter era: H(φ)=HF[1+(σ φ)^{-3}]^{1/2} (from eq.(62) related forms)
    // Actually from text: H(φ) = HF / [1+(σ φ)^{-3}]^{1/2} if rearranged.

    let sigma_f = Float::with_val(128, sigma);
    let one_minus_nu = Float::with_val(128, 1.0 - NU);
    let hf_f = Float::with_val(128, HF);
    let phi_val = phi.clone();
    let phidot_val = phidot.clone();

    let y = (sigma_f.clone()*&phi_val).pow(-3);
    let h_val = &hf_f*(Float::with_val(128,1.0)+&y).sqrt();

    let phidot_new = &one_minus_nu*&phi_val*&h_val;

    // Derivatives:
    // dH/dφ = hf d/dφ[(1+y)^{1/2}] with y=(σφ)^{-3}
    // dy/dφ = -3 σ^{-3} φ^{-4}, careful:
    // Actually y=(σ φ)^{-3}, dy/dφ = -3 (σ^(-3)) φ^(-4)? No σ^(-3) not correct since σ>0:
    // y=(σφ)^{-3}=σ^{-3}φ^{-3}, dy/dφ=-3σ^{-3}φ^{-4}
    let dy_dphi = Float::with_val(128, -3.0)*sigma_f.pow(-3)*phi_val.pow(-4);
    let dH_dphi = &hf_f*Float::with_val(128,0.5)*(Float::with_val(128,1.0)+&y).pow(-1.0/2.0)*dy_dphi;

    let dphidot_dphi = &one_minus_nu*( &h_val + &phi_val*&dH_dphi );
    let phiddot = &dphidot_dphi*&phidot_new;

    let phi_next = phi_val + &phidot_val*Float::with_val(128, DT);
    let phidot_next = phidot_val + phiddot*Float::with_val(128, DT);

    (phi_next, phidot_next)
}

fn main() {
    // Initialize φ and φ̇ at early times. The universe starts at primeval de Sitter vacuum: H=HI, φ ~ 0.
    // At φ=0, H=HI, we can choose φ small and phidot based on eq.(56) for radiation era start.
    let alpha = alpha_param();
    let sigma = sigma_param();

    let mut phi = Float::with_val(128, 1e-30); // very small initial φ
    let h_initial = HI/(1.0_f64+(alpha*1e-30).powi(4)).sqrt();
    let phidot_init = (1.0 - NU)*1e-30*h_initial;
    let mut phidot = Float::with_val(128, phidot_init);

    // We evolve from early times (radiation) to matter era and towards final de Sitter.
    // Transition criteria: once we approach conditions for matter domination:
    // After radiation era: at φ=α^{-1}, end of inflation and radiation = matter equivalence.
    let phi_end_radiation = Float::with_val(128, 1.0/alpha);

    let mut file = File::create("quantum_gravity_plasma_data.csv").unwrap();
    writeln!(file, "step,phi,phidot,H,Potential,Kinetic,Era").unwrap();

    // Radiation-era evolution
    for step in 0..TIME_STEPS {
        let phi_f = phi.to_f64().unwrap();
        let pot = potential_radiation(phi_f, alpha);
        let kin = kinetic_radiation(phi_f, alpha);
        let x = (alpha*phi_f).powi(4);
        let h_val = HI/(1.0+x).sqrt();
        writeln!(file, "{},{},{},{},{},{},{}",
                 step, phi_f, phidot.to_f64().unwrap(),
                 h_val, pot, kin, "radiation").unwrap();

        if phi > phi_end_radiation {
            break;
        }
        let (phi_new, phidot_new) = evolve_radiation(&phi, &phidot, alpha);
        phi = phi_new;
        phidot = phidot_new;
    }

    // Now evolve matter era
    // For matter era: we continue from last φ and phidot.
    for step in TIME_STEPS..(2*TIME_STEPS) {
        let phi_f = phi.to_f64().unwrap();
        let pot = potential_matter(phi_f, sigma);
        let kin = kinetic_matter(phi_f, sigma);
        let y = (sigma*phi_f).powf(-3.0);
        let h_val = HF*(1.0+y).sqrt();

        writeln!(file, "{},{},{},{},{},{},{}",
                 step, phi_f, phidot.to_f64().unwrap(),
                 h_val, pot, kin, "matter").unwrap();

        let (phi_new, phidot_new) = evolve_matter(&phi, &phidot, sigma);
        phi = phi_new;
        phidot = phidot_new;

        // As we evolve, H→HF as φ→∞. We can stop after some large φ ensures final de Sitter.
        if h_val - HF.abs() < 1e-35 {
            break;
        }
    }

    println!("Simulation complete. Data in quantum_gravity_plasma_data.csv");
}
