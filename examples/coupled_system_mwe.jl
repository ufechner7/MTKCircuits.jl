# Minimal working example that reproduces:
# ┌ Warning: Did not converge after `maxiters = 100` substitutions.
# └ @ Symbolics ~/.julia/packages/Symbolics/.../src/variable.jl:451
#
# Key elements required:
#   - Single-phase PMSG with algebraic torque variable Tₑ
#   - Nested DiodeBridge subsystem (2 diodes + 2 bleed resistors)
#   - Capacitor across the DC load
#   - Two leak resistors (one on each DC rail to Ground)

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

# Single Shockley diode
function ShockleyDiode(; name, I_S = 1.0e-6, n_ideal = 1.0, T_K = 293.15)
    @parameters I_s=I_S n_id=n_ideal T=T_K k_B=1.380649e-23 q_e=1.602e-19
    @variables v(t)=0.0 i(t)
    systems = @named begin p = Pin(); n = Pin() end
    eqs = [
        v ~ p.v - n.v,
        i ~ I_s * (exp(clamp(v / (n_id * k_B * T / q_e), -40, 40)) - 1),
        p.i ~ i,
        n.i ~ -i,
    ]
    return System(eqs, t; name, systems)
end

# Single-phase full-wave bridge: 2 diodes + 2 bleed resistors (nested subsystem)
function DiodeBridge(; name)
    @named D1 = ShockleyDiode()
    @named D2 = ShockleyDiode()
    @named R1 = Resistor(R = 1.0e8)
    @named R2 = Resistor(R = 1.0e8)
    eqs = [
        connect(D1.p, D2.n),
        connect(D1.p, R1.p), connect(D1.n, R1.n),
        connect(D2.p, R2.p), connect(D2.n, R2.n),
    ]
    return System(eqs, t; systems = [D1, D2, R1, R2], name)
end

# Single-phase permanent-magnet synchronous generator
function PMSG(; name, R = 0.05, L = 0.000635, Ψ = 0.192, J = 0.011, F = 0.001889, p = 4, Tf = 0.0)
    @named input1 = RealInput()
    @named pin_a = Pin()
    @named pin_n = Pin()
    @parameters R=R L=L Ψ=Ψ J=J F=F p=p Tf=Tf
    @variables ia(t)=0 ωₘ(t)=0 Θ(t)=0 Tₑ(t)=0
    D = Differential(t)
    eqs = [
        D(ia)  ~ (1/L) * (Ψ*p*ωₘ*sin(p*Θ) - R*ia - (pin_a.v - pin_n.v)),
        Tₑ     ~ Ψ*p*ia*sin(p*Θ),
        D(ωₘ)  ~ (1/J) * (input1.u - Tₑ - Tf - F*ωₘ),
        D(Θ)   ~ ωₘ,
        pin_a.i ~ -ia,
        pin_a.i + pin_n.i ~ 0,
    ]
    return System(eqs, t, [ia, ωₘ, Θ, Tₑ], [R, L, Ψ, J, F, p, Tf];
        systems = [input1, pin_a, pin_n], name)
end

@named aero   = Constant(k = -97)
@named gen    = PMSG(p = 4)
@named bridge = DiodeBridge()
@named rload  = Resistor(R = 10.0)
@named gnd    = Ground()
@named leak_p = Resistor(R = 1.0e9)
@named leak_n = Resistor(R = 1.0e9)
@named cap    = Capacitor(C = 1.0e-6, v = 0.0)

eqs = [
    connect(aero.output, gen.input1),
    connect(gen.pin_a, bridge.D1.p, bridge.D2.n),
    connect(rload.p, bridge.D1.n),
    connect(rload.n, bridge.D2.p),
    connect(rload.p, leak_p.p), connect(leak_p.n, gnd.g),
    connect(rload.n, leak_n.p), connect(leak_n.n, gnd.g),
    connect(cap.p, rload.p),    connect(cap.n, rload.n),
    connect(gen.pin_n, gnd.g),
]

@named sys = System(eqs, t; systems = [aero, gen, rload, leak_p, leak_n, bridge, cap, gnd])
sysc = mtkcompile(sys; warn_initialize_determined = false)
prob = ODEProblem(sysc, [], (0.0, 1.0); warn_initialize_determined = false)
nothing