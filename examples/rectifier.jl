using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using OrdinaryDiffEq
using Plots

begin
	function ShockleyDiode(; name, I_S = 1e-6, n_ideal = 1.0, T_K = 293.15)
	
	    @parameters begin
	        I_s   = I_S
	        n_id  = n_ideal
	        T     = T_K
	        k_B   = 1.380649e-23   
	        q_e   = 1.602e-19     
	    end
	    @variables begin
	        v(t) =0.0
	        i(t) 
	    end
	    systems = @named begin
	        p = Pin()
	        n = Pin()
	    end
	
	    eqs = [
	        v   ~ p.v - n.v,
	        #i   ~ I_s * (exp(v / (n_id * k_B * T / q_e)) - 1),
			i ~ I_s*(exp(clamp(v/(n_id*k_B*T/q_e),-40,40))-1),
	        p.i ~  i,
	        n.i ~ -i,
	    ]
	    System(eqs, t; name, systems)
	end

	

	
	function SourceAC3(; name,
	                          RMSVoltage = 220.0,
	                          Phase      = 0.0,
	                          Freq       = 50.0,
	                          τ          = 1e-6)
	    @named pin_a = Pin()
	    @named pin_b = Pin()
	    @named pin_c = Pin()
	    @named pin_n = Pin()
	
	    @parameters RMSVoltage=RMSVoltage  Phase=Phase  Freq=Freq  τ=τ
	
Vpeak = RMSVoltage * sqrt(2.0)
	    ramp  = 1 - exp(-t / τ)
	
	    eqs = [
	        pin_a.v - pin_n.v ~ ramp * Vpeak * sin(2π * Freq * t + Phase),
	        pin_b.v - pin_n.v ~ ramp * Vpeak * sin(2π * Freq * t + Phase - 2π/3),
	        pin_c.v - pin_n.v ~ ramp * Vpeak * sin(2π * Freq * t + Phase - 4π/3),
	
	        pin_a.i + pin_b.i + pin_c.i + pin_n.i ~ 0,
	    ]
	
	    System(eqs, t, [], [RMSVoltage, Phase, Freq, τ];
	        systems = [pin_a, pin_b, pin_c, pin_n],
	        name)
	end

	
	function DiodeBridge(; name, I_S = 1e-6, n_ideal = 1.0, T_K = 293.15)
	
	    @named D1    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
	    @named D2    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
		@named D3    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
		@named D4    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
		@named D5    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
		@named D6    = ShockleyDiode(I_S = I_S, n_ideal = n_ideal, T_K = T_K)
		@named R1 = Resistor(R=1e8)
		@named R2 = Resistor(R=1e8)
		@named R3 = Resistor(R=1e8)
		@named R4 = Resistor(R=1e8)
		@named R5 = Resistor(R=1e8)
		@named R6 = Resistor(R=1e8)
	
		eqs = [
		   
	    connect(D1.p, D2.n),
	    connect(D3.p, D4.n),
	    connect(D5.p, D6.n),		
	    connect(D1.n, D3.n, D5.n),
	    connect(D2.p, D4.p, D6.p),
		connect(D1.p,R1.p),
		connect(D1.n,R1.n),
		connect(D2.p,R2.p),
		connect(D2.n,R2.n),
		connect(D3.p,R3.p),
		connect(D3.n,R3.n),
		connect(D4.p,R4.p),
		connect(D4.n,R4.n),
		connect(D5.p,R5.p),
		connect(D5.n,R5.n),
		connect(D6.p,R6.p),
		connect(D6.n,R6.n),
	
		]
	
	    System(eqs, t; systems = [D1, D2, D3, D4,D5 ,D6, R1, R2, R3, R4 , R5, R6], name)
	end

	function CircuitBridge(; name,
                               RMSVoltage = 9.0,
                               Freq       = 50.0,
                               R_load     = 10.0,
                               I_S        = 1e-6,
                               n_ideal    = 1.0,
                               T_K        = 293.15)

    
    @named src = SourceAC3(RMSVoltage = RMSVoltage, Freq = Freq, Phase = 0.0)
	@named bridge = DiodeBridge()
	@named gnd   = Ground()

    @named rload = Resistor(R = R_load)


	eqs = [
	   
    connect(src.pin_a, bridge.D1.p, bridge.D2.n),
    connect(src.pin_b, bridge.D3.p, bridge.D4.n),
    connect(src.pin_c, bridge.D5.p, bridge.D6.n),		
    connect(rload.p, bridge.D1.n, bridge.D3.n, bridge.D5.n),
    connect(rload.n, bridge.D2.p, bridge.D4.p, bridge.D6.p),
    connect(src.pin_n, gnd.g),

	]

    System(eqs, t; systems = [src, bridge ,gnd, rload], name)
	end

	@named circuitbridge =  CircuitBridge(RMSVoltage = 9.0, Freq = 50.0, R_load = 100.0)
	
	sysBridge  = mtkcompile(circuitbridge)
	probBridge = ODEProblem(sysBridge, [], (0.0, 0.1))
	
	solBridge  = solve(probBridge,Rodas5P(),abstol=1e-9, reltol=1e-12)
	plot(solBridge, idxs = [circuitbridge.rload.p.i])
end