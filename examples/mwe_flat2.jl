# Flat (component-free) equivalent of coupled_system_mwe.jl
#
# Writing the equations directly avoids the nested component / pin-alias
# substitution chains that trigger:
#   ┌ Warning: Did not converge after `maxiters = 100` substitutions.
#   └ @ Symbolics …/src/variable.jl:451
#
# Circuit topology — single-phase full-bridge rectifier
# ────────────────────────────────────────────────────────
#   Generator (PMSG): AC terminals va (pin_a) and vb (pin_n), both floating
#   D1  : Shockley diode,  anode = va,    cathode = vp  (DC+ rail)
#   D2  : Shockley diode,  anode = gnd=0, cathode = va
#   D3  : Shockley diode,  anode = vb,    cathode = vp  (DC+ rail)
#   D4  : Shockley diode,  anode = gnd=0, cathode = vb
#   R1–R4: bleed resistors parallel to each diode
#   rload : DC load          (vp → gnd)
#   cap   : DC-bus capacitor (vp → gnd) — state v_cap = vp  (vn = gnd = 0)
#   leak_p: rail leakage     (vp → gnd)
#
# State variables    : ia, ωₘ, Θ, v_cap
# Algebraic variables: va, vb  (vp = v_cap, vn = 0)

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using ControlPlots

D = Differential(t)

# ── Parameters (plain Julia constants — no MTK substitution chains) ─────────────
const _R_g     = 0.05           # winding resistance           [Ω]
const _L_g     = 0.000635       # winding inductance            [H]
const _Ψ_g     = 0.192          # flux linkage                  [Wb]
const _J_g     = 0.011          # rotor inertia                 [kg·m²]
const _F_g     = 0.001889       # viscous damping               [N·m·s]
const _p_g     = 4              # pole pairs
const _T_aero  = -97.0/100          # aerodynamic torque (Constant k = -97)

const _I_s     = 1.0e-6                             # diode sat. current   [A]
const _V_T     = 1.380649e-23 * 293.15 / 1.602e-19  # thermal voltage ≈ 0.02527 V

const _R_load  = 10.0           # DC load resistance            [Ω]
const _R_bleed = 1.0e8          # diode bleed resistor          [Ω]
const _R_leak  = 1.0e9          # DC rail leakage to ground     [Ω]
const _C_cap   = 1.0e-6         # DC-bus capacitance            [F]

# ── Variables ────────────────────────────────────────────────────────────────────
@variables ia(t)=0.0  ωₘ(t)=0.0  Θ(t)=0.0  v_cap(t)=0.0   # differential states
@variables va(t)=0.0  vb(t)=0.0                               # algebraic node voltages (guesses)

# ── Symbolic branch-current expressions (inlined — no extra variables) ──────────
#   D1 (va → v_cap), D2 (0 → va), D3 (vb → v_cap), D4 (0 → vb);  vn = 0
_iD1 = _I_s * (exp(clamp((va - v_cap) / _V_T, -40.0, 40.0)) - 1)
_iD2 = _I_s * (exp(clamp(-va          / _V_T, -40.0, 40.0)) - 1)
_iD3 = _I_s * (exp(clamp((vb - v_cap) / _V_T, -40.0, 40.0)) - 1)
_iD4 = _I_s * (exp(clamp(-vb          / _V_T, -40.0, 40.0)) - 1)
_iR1 = (va - v_cap) / _R_bleed
_iR2 = -va          / _R_bleed
_iR3 = (vb - v_cap) / _R_bleed
_iR4 = -vb          / _R_bleed

eqs = [
    # ── PMSG dynamics ──────────────────────────────────────────────────────────
    # Armature loop voltage: coil driven by (va − vb), both terminals floating
    D(ia)    ~ (1 / _L_g) * (_Ψ_g * _p_g * ωₘ * sin(_p_g * Θ) - _R_g * ia - (va - vb)),
    # Newton's law for rotor (Tₑ = Ψ·p·ia·sin(p·Θ) inlined)
    D(ωₘ)   ~ (1 / _J_g) * (_T_aero - _Ψ_g * _p_g * ia * sin(_p_g * Θ) - _F_g * ωₘ),
    D(Θ)    ~ ωₘ,

    # ── Capacitor ODE: D(v_cap) = i_cap / C ────────────────────────────────────
    # i_cap = (iD1+iR1) + (iD3+iR3) − v_cap/R_load − v_cap/R_leak
    D(v_cap) ~ (_iD1 + _iR1 + _iD3 + _iR3 - v_cap / _R_load - v_cap / _R_leak) / _C_cap,

    # ── KCL at va (pin_a) →  determines va algebraically ───────────────────────
    # ia leaves va via D1/R1; iD2/iR2 enter va from gnd
    ia - _iD1 - _iR1 + _iD2 + _iR2 ~ 0,

    # ── KCL at vb (pin_n, floating) →  determines vb algebraically ─────────────
    # ia enters vb (return path); iD3/iR3 leave vb to vp; iD4/iR4 enter vb from gnd
    -ia - _iD3 - _iR3 + _iD4 + _iR4 ~ 0,
]

@named sys = System(eqs, t)
sysc = mtkcompile(sys; warn_initialize_determined = false)
prob = ODEProblem(sysc, [], (0.0, 3.0); warn_initialize_determined = false)

sol = solve(prob, Rodas5P(), abstol = 1e-9, reltol = 1e-12, dense = false, saveat = 0.0001, maxiters = 10000000)

time    = sol.t
ω_rpm   = sol[ωₘ] .* (30 / π)           # rad/s → RPM
v_dc    = sol[v_cap]                      # DC-bus voltage (= vp − vn)
i_a     = sol[ia]                         # armature current
i_load  = v_dc ./ _R_load                 # DC load current

plotx(time, ω_rpm, v_dc, i_load, i_a;
      ylabels = ["Speed [RPM]", "DC voltage [V]", "Load current [A]", "Armature current [A]"],
      labels  = ["ωₘ", "v_dc", "i_load", "ia"])
