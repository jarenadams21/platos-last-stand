use std::fs::File;
use std::io::Write;
use std::f64::consts::PI;

// We assume a monoatomic ideal gas with degrees of freedom f=3, γ = Cp/Cv = 5/3 for demonstration.
// For generality, we keep constants dimensionless. k_B = 1, N = 1. Thus, P V = T.

// Carnot cycle parameters
const T_H: f64 = 500.0;   // Hot reservoir temperature (K), example
const T_C: f64 = 300.0;   // Cold reservoir temperature (K), example
const GAMMA: f64 = 5.0/3.0; // Ratio of specific heats for monoatomic ideal gas
const N_STEPS_ISOTHERMAL: usize = 200; // Number of steps for isothermal processes
const N_STEPS_ADIABATIC: usize = 200;  // Number of steps for adiabatic processes

// Initial volumes and conditions
// Let's choose an initial state A:
// We start at (P1, V1, T_H). For convenience, pick V1 = 1.0 arbitrary units.
// Then from T_H and P V = T, we get P1 = T_H / V1 = T_H.
const V1: f64 = 1.0;
const T_START: f64 = T_H;

// For the Carnot cycle, we must choose a volume ratio so that after the isothermal expansion at T_H,
// the gas expands from V1 to V2. We can pick V2 so that the adiabatic brings it down to T_C.
// Adiabatic relation: T V^(γ-1) = const. Since T_H * V1^(γ-1) = T_C * V2^(γ-1) for the second step reverse,
// We can solve for V2:
// T_H * V1^(γ-1) = T_C * V2^(γ-1) => V2 = V1 * (T_H/T_C)^(1/(γ-1))
fn v2() -> f64 {
    V1 * (T_H/T_C).powf(1.0/(GAMMA-1.0))
}
// After adiabatic expansion to state C at T_C, we do isothermal compression at T_C back down to volume V3.
// Efficiency requires a certain ratio. For a Carnot cycle, V3 = V2 * (T_C/T_H)^(1/(γ-1)) would bring back to T_C after adiabatic compression.
// But we know at the end we must return to the original state A.
// Actually, by symmetry of the Carnot cycle in T-S space:
// - After the isothermal at T_C (C->D), we must have the same ratio of volumes as the hot isothermal, so that the final adiabatic (D->A) returns to T_H.
// Let V2_val = v2():
// On the isothermal at T_C, we choose a compression ratio so that after the final adiabatic we return to V1 at T_H.
// The final V4 must be equal to the initial V1.
// From D->A adiabatic: T_C * V3^(γ-1) = T_H * V4^(γ-1).
// We know V4 = V1. Thus, V3 = V1 * (T_H/T_C)^(1/(γ-1)) which is exactly symmetric.
// Therefore, V3 = V2 and V4 = V1 for a perfect Carnot cycle! 
// Actually, for a Carnot cycle: V2 and V3 differ. After isothermal expansion at T_H, volume is larger.
// After adiabatic expansion, volume is even larger.
// Let's re-derive carefully:

// Carnot cycle states:
// A(T_H,V1) -> Isothermal expansion at T_H to B(T_H,V2).
// Adiabatic expansion: B(T_H,V2) -> C(T_C,V3).
// Isothermal compression at T_C: C(T_C,V3) -> D(T_C,V4).
// Adiabatic compression: D(T_C,V4) -> A(T_H,V1).
//
// Adiabatic relations:
// T_H * V2^(γ-1) = T_C * V3^(γ-1) => V3 = V2 * (T_H/T_C)^(1/(γ-1)).
// Similarly, for the last adiabatic:
// T_C * V4^(γ-1) = T_H * V1^(γ-1) => V4 = V1 * (T_H/T_C)^(1/(γ-1)).
//
// We started with V2 = V1 * (T_H/T_C)^(1/(γ-1)) from the first adiabatic relation inverted. Let's correct that:
// Actually, for the first isothermal at T_H, we just pick a ratio. Let's define a compression ratio r:
// Let’s choose a convenient ratio for isothermal expansion: say we double the volume from V1=1.0 to V2=2.0.
// Then from B to C (adiabatic): T_H * V2^(γ-1) = T_C * V3^(γ-1).
// Solve for V3: V3 = V2 * (T_H/T_C)^(1/(γ-1)).
// Next, from C to D (isothermal at T_C), we must compress. Let's say we compress until some volume V4 < V3. 
// The final adiabatic must bring us back to T_H and V1=1.0:
// T_C * V4^(γ-1) = T_H * V1^(γ-1) => V4 = V1 * (T_H/T_C)^(1/(γ-1)) = (T_H/T_C)^(1/(γ-1)).
  
// For a perfect Carnot cycle symmetry, let's set it so that the adiabatic expansions and compressions are symmetric. 
// We want a neat, closed loop. Let's pick initial conditions to simplify the math:

// Let's pick V2 so that after the adiabatic to T_C we have a nice ratio:
// Let’s define V2 = V1 * 2.0 (just a choice).
// Then from T_H to T_C adiabatic: 
// V3 = V2 * (T_H/T_C)^(1/(γ-1)).

// For the T_C isothermal, we compress from V3 down to V4. After final adiabatic, we must return to V1=1.0:
// So from D->A adiabatic: V4^(γ-1) = (T_H/T_C)*V1^(γ-1) = (T_H/T_C)*1^(γ-1) = T_H/T_C
// V4 = (T_H/T_C)^(1/(γ-1)).

// To close the loop, we must have V4 < V1 * 2.0 to ensure a proper cycle. Let's verify numerical values:

// Let’s compute with T_H=500, T_C=300, γ=5/3:
// (γ-1)=2/3
// (T_H/T_C)^(1/(γ-1)) = (500/300)^(3/2) = (5/3)^(1.5) ~ (1.6667)^1.5 ≈ 2.15 approx.
// If V2=2.0, then V3 = V2 * (T_H/T_C)^(1/(γ-1)) = 2.0 * 2.15 = 4.3 approx.
// V4 = (T_H/T_C)^(1/(γ-1)) = 2.15 approx.

// So the cycle is: A(500K,1.0), B(500K,2.0), C(300K,4.3), D(300K,2.15), A(500K,1.0).
// This forms a closed loop. Perfectly fine for demonstration.

// We'll do stepwise transformations:
// - From A to B: isothermal at T_H, we go from V=1.0 to V=2.0 in N_STEPS_ISOTHERMAL steps
// - From B to C: adiabatic, from V=2.0 to V=4.3 in N_STEPS_ADIABATIC steps
// - From C to D: isothermal at T_C, from V=4.3 to V=2.15 in N_STEPS_ISOTHERMAL steps
// - From D to A: adiabatic, from V=2.15 to V=1.0 in N_STEPS_ADIABATIC steps

// We'll record at each step: step_number, V, P, T, S, Q (cumulative), W (cumulative), "phase".
// Entropy S for ideal gas (up to constants): S = const + ln(V) + (f/2)*ln(T), 
// but since we only need changes, we can set S0 at A as reference. 
// Let's just track S changes. For a monoatomic ideal gas: S ∝ ln(V T^{3/2}) (since f=3).
// S difference: ΔS = ln(V) + (3/2)*ln(T) - [ln(V_A)+ (3/2)*ln(T_A)]
// We'll just store relative S with respect to the initial point A where S(A)=0.

// Heat Q is absorbed (Q>0) if dQ>0 on isothermal expansions and negative on isothermal compressions.
// On adiabatic steps, Q=0 change. Work W = ∫P dV. On isothermal steps, Q = W. On adiabatic steps, W changes U and no Q exchange.

// Efficiency check: η = 1 - T_C/T_H. With T_H=500, T_C=300 => η=1-0.6=0.4=40%.

fn main() {
    let mut file = File::create("carnot_data.csv").unwrap();
    writeln!(file, "step,V,P,T,S,Q,W,phase").unwrap();

    let v2 = 2.0;
    let ratio = (T_H/T_C).powf(1.0/(GAMMA-1.0));
    let v3 = v2 * ratio;
    let v4 = ratio;

    // Initial State A:
    let mut v = V1;
    let mut t = T_H;
    let p = t/v;
    // Reference entropy at A = 0
    let s_ref = (v).ln() + 1.5*(t).ln(); 
    let mut s = 0.0;
    let mut q_cumulative = 0.0;
    let mut w_cumulative = 0.0;
    let mut step_count = 0;

    // Helper to write data
    let mut write_data = |step: usize, v_val: f64, t_val: f64, q_val: f64, w_val: f64, phase: &str, file: &mut File| {
        let p_val = t_val / v_val;
        let s_val = (v_val).ln() + 1.5*(t_val).ln() - s_ref;
        writeln!(file, "{},{},{},{},{},{},{},{}",
                 step, v_val, p_val, t_val, s_val, q_val, w_val, phase).unwrap();
    };

    // 1) Isothermal expansion at T_H: from V1 to V2
    // dW = P dV = (T_H/V) dV, dQ = dW since ΔU=0 in isothermal.
    // Integrate in small steps:
    let dv = (v2 - v)/(N_STEPS_ISOTHERMAL as f64);
    for i in 0..N_STEPS_ISOTHERMAL {
        let p_local = t/v;
        let dW = p_local * dv;
        let dQ = dW;
        w_cumulative += dW;
        q_cumulative += dQ;
        v += dv; // increment volume
        step_count += 1;
        write_data(step_count, v, t, q_cumulative, w_cumulative, "isothermal_hot", &mut file);
    }

    // Now at point B(T_H,V2), do adiabatic expansion B->C:
    // Adiabatic: T V^(γ-1) = const
    // Let's find the constant: T * V^(γ-1) = const
    let c_adiab1 = t * v.powf(GAMMA-1.0);

    let dv = (v3 - v)/(N_STEPS_ADIABATIC as f64);
    for _ in 0..N_STEPS_ADIABATIC {
        v += dv;
        // T = c_adiab1 / V^(γ-1)
        t = c_adiab1 / v.powf(GAMMA-1.0);
        // dW = ∫P dV, P = T/V
        // But we must do stepwise: For a small step dW = P * dV
        // Q=0 here, so all energy change is work changing internal energy
        let p_local = t/v;
        let dW = p_local * dv;
        w_cumulative += dW;
        // Q=0 in adiabatic
        step_count += 1;
        write_data(step_count, v, t, q_cumulative, w_cumulative, "adiabatic_expand", &mut file);
    }

    // At C(T_C,V3), isothermal compression at T_C: from V3 to V4
    // T = T_C constant
    t = T_C;
    let dv = (v4 - v)/(N_STEPS_ISOTHERMAL as f64);
    for _ in 0..N_STEPS_ISOTHERMAL {
        v += dv;
        let p_local = t/v;
        let dW = p_local * dv; 
        // Since this is compression, dv < 0 actually. Wait, we must compress from larger volume to smaller volume,
        // so v4 < v3? We must check if v4 < v3:
        // v3 ~4.3, v4 ~2.15, indeed v4 < v3, so dv is negative in code if we didn't reorder.
        // Let's correct: v4 < v3, we must go from v3 to v4, so dv should be negative if we do v4 - v3.
        // Let's redefine dv after we know v3 and v4 properly.

        // Correcting the sign and approach:
    }

    // We made a mistake: we defined dv = (v4 - v)/(N_STEPS_ISOTHERMAL)
    // Since v4 < v3, this will be negative. That's actually correct because we are compressing.
    // When we add dv to v, v will decrease. Just keep careful track of sign:
    let mut v_current = v3;
    let dv_iso_cold = (v4 - v3)/(N_STEPS_ISOTHERMAL as f64);
    for _ in 0..N_STEPS_ISOTHERMAL {
        let p_local = T_C / v_current;
        // dW = p dV. Note dv is negative. W_cumulative changes by p*dV (dV<0)
        let dW = p_local * dv_iso_cold;
        // On isothermal compression at T_C, Q = -W (since ΔU=0).
        // If we do positive from system perspective:
        // Actually, when we compress, we do work on the gas, W>0 from environment. 
        // The gas loses heat Q = W (still from system perspective, if system's W is positive out, now it's reversed).
        // Careful with sign conventions:
        // In thermodynamics: W>0 if work done by system. Here, compression is done on the system, so from system perspective W<0.
        // dV<0 => p*dV<0 => dW<0 from system perspective (system is not doing work, work is done on it).
        // Q = W for isothermal, so Q<0 as well, system releases heat to cold reservoir.
        w_cumulative += dW;
        q_cumulative += dW; // Q = W in isothermal

        v_current += dv_iso_cold;
        step_count += 1;
        write_data(step_count, v_current, T_C, q_cumulative, w_cumulative, "isothermal_cold", &mut file);
    }

    // Update v, t at the end of this step
    v = v_current;
    t = T_C;

    // Adiabatic compression D->A:
    // We know final state: V1=1.0, T_H=500K
    // Adiabatic relation again: T * V^(γ-1) = const
    let c_adiab2 = t * v.powf(GAMMA-1.0);
    let dv_adiab_back = (V1 - v)/(N_STEPS_ADIABATIC as f64);
    for _ in 0..N_STEPS_ADIABATIC {
        v += dv_adiab_back;
        t = c_adiab2 / v.powf(GAMMA-1.0);
        let p_local = t/v;
        let dW = p_local * dv_adiab_back;
        w_cumulative += dW;
        // Q=0 adiabatic
        step_count += 1;
        write_data(step_count, v, t, q_cumulative, w_cumulative, "adiabatic_compress", &mut file);
    }

    println!("Simulation complete. Data in carnot_data.csv");
}