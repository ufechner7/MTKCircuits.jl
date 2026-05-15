using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots

begin
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

    function LoadABC(; name, R_L = 10.0)

        @named pin_a = Pin()
        @named pin_b = Pin()
        @named pin_c = Pin()
        @named pin_n = Pin()

        @parameters R_L = R_L

        eqs = [

            pin_a.v - pin_n.v ~ R_L * pin_a.i,
            pin_b.v - pin_n.v ~ R_L * pin_b.i,
            pin_c.v - pin_n.v ~ R_L * pin_c.i,


            pin_n.i ~ -(pin_a.i + pin_b.i + pin_c.i),
        ]

        return System(
            eqs, t, [], [R_L];
            systems = [pin_a, pin_b, pin_c, pin_n],
            name
        )
    end

    function WindResistive(; name, R_L = 10.0, p = 4)

        @named aero = Constant(k = -97.24)

        @named gen = PMSG(p = p)
        @named Aload = LoadABC(R_L = R_L)
        @named gnd = Ground()

        eqs = [

            connect(aero.output, gen.input1),


            connect(gen.pin_a, Aload.pin_a),
            connect(gen.pin_b, Aload.pin_b),
            connect(gen.pin_c, Aload.pin_c),


            connect(gen.pin_n, gnd.g),
        ]

        return System(
            eqs, t;
            systems = [aero, gen, Aload, gnd],
            name
        )
    end

    @named SystemWindResistive = WindResistive(R_L = 10.0)

    sysWindResistive = mtkcompile(SystemWindResistive)

    probWindResistive = ODEProblem(sysWindResistive, [], (0.0, 15.0))
    solWindResistive = solve(probWindResistive, Rodas5P(), dt = 0.00001, adaptive = false, maxiters = Int(1.0e9))
    plot(solWindResistive, idxs = [SystemWindResistive.Aload.pin_a.i], tspan = (0, 10))
end
