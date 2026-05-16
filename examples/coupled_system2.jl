using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using OrdinaryDiffEqFIRK: RadauIIA5
using ControlPlots

function ShockleyDiode(; name, I_S = 1.0e-6, n_ideal = 1.0, T_K = 293.15)

    @parameters begin
        I_s = I_S
        n_id = n_ideal
        T = T_K
        k_B = 1.380649e-23
        q_e = 1.602e-19
    end
    @variables begin
        v(t) = 0.0
        i(t)
    end
    
    systems = @named begin
        p = Pin()
        n = Pin()
    end

    eqs = [
        v ~ p.v - n.v,
        #i   ~ I_s * (exp(v / (n_id * k_B * T / q_e)) - 1),
        i ~ I_s * (exp(clamp(v / (n_id * k_B * T / q_e), -40, 40)) - 1),
        p.i ~ i,
        n.i ~ -i,
    ]
    return System(eqs, t; name, systems)
end

function DiodeBridge(; name, I_S = 1.0e-6, n_ideal = 1.0, T_K = 293.15)

    @named D1 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named D2 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named D3 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named D4 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named D5 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named D6 = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
    @named C1 = Capacitor(C = 625e-12, v = 0.0)
    @named C2 = Capacitor(C = 625e-12, v = 0.0)
    @named C3 = Capacitor(C = 625e-12, v = 0.0)
    @named C4 = Capacitor(C = 625e-12, v = 0.0)
    @named C5 = Capacitor(C = 625e-12, v = 0.0)
    @named C6 = Capacitor(C = 625e-12, v = 0.0)
    @named R1 = Resistor(R = 1.0e8)
    @named R2 = Resistor(R = 1.0e8)
    @named R3 = Resistor(R = 1.0e8)
    @named R4 = Resistor(R = 1.0e8)
    @named R5 = Resistor(R = 1.0e8)
    @named R6 = Resistor(R = 1.0e8)
    @named RS1 = Resistor(R = 0.001)
    @named RS2 = Resistor(R = 0.001)
    @named RS3 = Resistor(R = 0.001)
    @named RS4 = Resistor(R = 0.001)
    @named RS5 = Resistor(R = 0.001)
    @named RS6 = Resistor(R = 0.001)

    eqs = [

        connect(D1.p, D2.n),
        connect(D3.p, D4.n),
        connect(D5.p, D6.n),
        connect(D1.n, RS1.p),
        connect(D3.n, RS3.p),
        connect(D5.n, RS5.p),
        connect(RS1.n, RS3.n, RS5.n),
        connect(RS2.n, D2.p),
        connect(RS4.n, D4.p),
        connect(RS6.n, D6.p),
        connect(RS2.p, RS4.p, RS6.p),
        connect(D1.p, R1.p),
        connect(D1.n, R1.n),
        connect(D2.p, R2.p),
        connect(D2.n, R2.n),
        connect(D3.p, R3.p),
        connect(D3.n, R3.n),
        connect(D4.p, R4.p),
        connect(D4.n, R4.n),
        connect(D5.p, R5.p),
        connect(D5.n, R5.n),
        connect(D6.p, R6.p),
        connect(D6.n, R6.n),
        connect(D1.p, C1.p),
        connect(D1.n, C1.n),
        connect(D2.p, C2.p),
        connect(D2.n, C2.n),
        connect(D3.p, C3.p),
        connect(D3.n, C3.n),
        connect(D4.p, C4.p),
        connect(D4.n, C4.n),
        connect(D5.p, C5.p),
        connect(D5.n, C5.n),
        connect(D6.p, C6.p),
        connect(D6.n, C6.n),

    ]

    return System(eqs, t; systems = [D1, D2, D3, D4, D5, D6, C1, C2, C3, C4, C5, C6, R1, R2, R3, R4, R5, R6, RS1, RS2, RS3, RS4, RS5, RS6], name)
end

function PMSG(;
        name,
        R = 0.05,
        L = 0.000635,
        Ψ = 0.192,
        J = 0.011,
        F = 0.001889,
        p = 4,
        Tf = 0.0
    )

    @named input1 = RealInput()
    @named pin_a = Pin()
    @named pin_b = Pin()
    @named pin_c = Pin()
    @named pin_n = Pin()

    @named output1 = RealOutput()
    @named output2 = RealOutput()

    @parameters R = R L = L Ψ = Ψ J = J F = F p = p Tf = Tf
    @variables  ia(t) = 0 ib(t) = 0 ic(t) = 0 ωₘ(t) = 0 Θ(t) = 0 Tₑ(t) = 0
    D = Differential(t)

    eqs = [

        D(ia) ~ (1 / L) * (Ψ * p * ωₘ * sin(p * Θ) - R * ia - (pin_a.v - pin_n.v)),
        D(ib) ~ (1 / L) * (Ψ * p * ωₘ * sin(p * Θ - 2π / 3) - R * ib - (pin_b.v - pin_n.v)),
        D(ic) ~ (1 / L) * (Ψ * p * ωₘ * sin(p * Θ + 2π / 3) - R * ic - (pin_c.v - pin_n.v)),

        Tₑ ~ Ψ * p * (ia * sin(p * Θ) + ib * sin(p * Θ - 2π / 3) + ic * sin(p * Θ + 2π / 3)),

        D(ωₘ) ~ (1 / J) * (input1.u - Tₑ - Tf - F * ωₘ),
        D(Θ) ~ ωₘ,

        pin_a.i ~ -ia,
        pin_b.i ~ -ib,
        pin_c.i ~ -ic,
        pin_a.i + pin_b.i + pin_c.i + pin_n.i ~ 0,

        output1.u ~ Tₑ,
        output2.u ~ ωₘ,
    ]

    return System(
        eqs, t,
        [ia, ib, ic, ωₘ, Θ, Tₑ],
        [R, L, Ψ, J, F, p, Tf];
        systems = [
            input1, pin_a, pin_b, pin_c, pin_n,
            output1, output2,
        ],
        name
    )
end

begin
    function WindDiodes(; name, R_L = 10.0, p = 4)

        @named aero = Constant(k = -97)
        #@named aero = Ramp(offset=0.0, height=-97.24, duration=15.0, start_time=0.0)
        @named gen = PMSG(p = p)
        @named zbridge = DiodeBridge()
        @named rload = Resistor(R = R_L)
        @named gnd = Ground()   # fija V_neutro = 0
        @named leak_p = Resistor(R = 1e6)
        @named leak_n = Resistor(R = 1e6)
        @named cap = Capacitor(C = 1e-6, v = 0.0)

        eqs = [

            connect(aero.output, gen.input1),

            connect(gen.pin_a, zbridge.D1.p, zbridge.D2.n),
            connect(gen.pin_b, zbridge.D3.p, zbridge.D4.n),
            connect(gen.pin_c, zbridge.D5.p, zbridge.D6.n),
            connect(rload.p, zbridge.RS1.n, zbridge.RS3.n, zbridge.RS5.n),
            connect(rload.n, zbridge.RS2.p, zbridge.RS4.p, zbridge.RS6.p),
            connect(rload.p, leak_p.p), connect(leak_p.n, gnd.g),
            connect(rload.n, leak_n.p), connect(leak_n.n, gnd.g),
            connect(cap.p, rload.p),
            connect(cap.n, rload.n),

            connect(gen.pin_n, gnd.g),
        ]

        return System(
            eqs, t;
            systems = [aero, gen, rload, leak_p, leak_n, zbridge, cap, gnd],
            name
        )
    end

end

@named sys = WindDiodes(R_L = 10.0, p = 4)

sysc = mtkcompile(sys; warn_initialize_determined = false)
prob = ODEProblem(sysc, [], (0.0, 1.0); warn_initialize_determined = false)

sol = solve(prob, RadauIIA5(κ = 0.005), abstol = 1e-7, reltol = 1e-8, saveat = 0.000025, maxiters = 5e6)

time    = sol.t
ω_rpm   = sol[sys.gen.ωₘ] .* (30 / π)    # rad/s → RPM
v_dc    = sol[sys.cap.v]                   # DC-bus voltage
i_load  = sol[sys.rload.i]                 # DC load current
i_a     = sol[sys.gen.ia]                  # phase-a armature current
i_b     = sol[sys.gen.ib]                  # phase-b armature current
i_c     = sol[sys.gen.ic]                  # phase-c armature current

plotx(time, ω_rpm, v_dc, i_load, i_a, i_b, i_c;
      ylabels = ["Speed [RPM]", "DC voltage [V]", "Load current [A]",
                 "ia [A]", "ib [A]", "ic [A]"],
      labels  = ["ωₘ", "v_dc", "i_load", "ia", "ib", "ic"])
