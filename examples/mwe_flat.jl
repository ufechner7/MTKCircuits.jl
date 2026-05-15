# Flat (component-free) equivalent of coupled_system_mwe.jl
#
# Writing the equations directly avoids the nested component / pin-alias
# substitution chains that trigger:
#   ┌ Warning: Did not converge after `maxiters = 100` substitutions.
#   └ @ Symbolics …/src/variable.jl:451
#
# Circuit topology
# ─────────────────
#   Generator (PMSG): AC terminal va, neutral at gnd (0 V)
#   D1  : Shockley diode,  anode = va,  cathode = vp  (DC+ rail)
#   D2  : Shockley diode,  anode = vn,  cathode = va  (DC− rail)
#   R1  : bleed resistor parallel to D1  (va → vp)
#   R2  : bleed resistor parallel to D2  (vn → va)
#   rload : DC load          (vp → vn)
#   cap   : DC-bus capacitor (vp → vn) — state v_cap = vp − vn
#   leak_p: rail leakage     (vp → gnd)
#   leak_n: rail leakage     (vn → gnd)
#
# State variables  : ia, ωₘ, Θ, v_cap
# Algebraic variables: va, vp   (vn = vp − v_cap is eliminated)

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq

D = Differential(t)

# ── Parameters (plain Julia constants — no MTK substitution chains) ─────────────
const _R_g     = 0.05           # winding resistance           [Ω]
const _L_g     = 0.000635       # winding inductance            [H]
const _Ψ_g     = 0.192          # flux linkage                  [Wb]
const _J_g     = 0.011          # rotor inertia                 [kg·m²]
const _F_g     = 0.001889       # viscous damping               [N·m·s]
const _p_g     = 4              # pole pairs
const _T_aero  = -97.0          # aerodynamic torque (Constant k = -97)

const _I_s     = 1.0e-6                             # diode sat. current   [A]
const _V_T     = 1.380649e-23 * 293.15 / 1.602e-19  # thermal voltage ≈ 0.02527 V

const _R_load  = 10.0           # DC load resistance            [Ω]
const _R_bleed = 1.0e8          # diode bleed resistor          [Ω]
const _R_leak  = 1.0e9          # DC rail leakage to ground     [Ω]
const _C_cap   = 1.0e-6         # DC-bus capacitance            [F]

# ── Variables ────────────────────────────────────────────────────────────────────
@variables ia(t)=0.0  ωₘ(t)=0.0  Θ(t)=0.0  v_cap(t)=0.0   # differential states
@variables va(t)  vp(t)                                       # algebraic node voltages

# ── Symbolic branch-current expressions (inlined — no extra variables) ──────────
#   D1 (va → vp),  D2 (vn → va),  R1 (va → vp),  R2 (vn → va);  vn = vp − v_cap
_iD1 = _I_s * (exp(clamp((va - vp)              / _V_T, -40.0, 40.0)) - 1)
_iD2 = _I_s * (exp(clamp((vp - v_cap - va)      / _V_T, -40.0, 40.0)) - 1)
_iR1 = (va - vp)         / _R_bleed
_iR2 = (vp - v_cap - va) / _R_bleed

eqs = [
    # ── PMSG dynamics ──────────────────────────────────────────────────────────
    # Faraday / Kirchhoff voltage law for armature loop (pin_n at gnd → va − 0 = va)
    D(ia)    ~ (1 / _L_g) * (_Ψ_g * _p_g * ωₘ * sin(_p_g * Θ) - _R_g * ia - va),
    # Newton's law for rotor (Tₑ = Ψ·p·ia·sin(p·Θ) inlined)
    D(ωₘ)   ~ (1 / _J_g) * (_T_aero - _Ψ_g * _p_g * ia * sin(_p_g * Θ) - _F_g * ωₘ),
    D(Θ)    ~ ωₘ,

    # ── Capacitor ODE: D(v_cap) = i_cap / C ────────────────────────────────────
    # i_cap comes from KCL at DC+ rail:  i_cap = iD1 + iR1 − v_cap/R_load − vp/R_leak
    D(v_cap) ~ (_iD1 + _iR1 - v_cap / _R_load - vp / _R_leak) / _C_cap,

    # ── KCL at AC node va  →  determines va algebraically ──────────────────────
    # ia (from generator) − iD1 + iD2 − iR1 + iR2 = 0
    ia - _iD1 + _iD2 - _iR1 + _iR2 ~ 0,

    # ── KCL at DC+ rail + DC− rail (sum)  →  determines vp algebraically ───────
    # iD1 + iR1 − vp/R_leak − iD2 − iR2 − (vp−v_cap)/R_leak = 0
    _iD1 + _iR1 - vp / _R_leak - _iD2 - _iR2 - (vp - v_cap) / _R_leak ~ 0,
]

@named sys = System(eqs, t)
sysc = mtkcompile(sys; warn_initialize_determined = false)
prob = ODEProblem(sysc, [], (0.0, 1.0); warn_initialize_determined = false)
nothing
