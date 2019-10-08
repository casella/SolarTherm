within SolarTherm.Models.CSP.CRS.Receivers;
model ParticleReceiver_1D
	extends Interfaces.Models.ReceiverFluid;
	
	parameter SI.Length H_rcv = 24.37 "Receiver drop height" annotation(Dialog(group="Technical data"));
	parameter SI.Length W_rcv = 24.37 "Receiver width" annotation(Dialog(group="Technical data"));
	parameter SI.Length L_rcv = 1 "Receiver length(depth)" annotation(Dialog(group="Technical data"));

	parameter SI.Area A = H_rcv * W_rcv "Receiver aperture area" annotation(Dialog(group="Technical data"));
	
	SI.Length th_c_in(start=0.43) "Curtain thicknesss at the inlet";

	parameter SI.Efficiency em = 0.86 "Emissivity" annotation(Dialog(group="Technical data"));
	parameter SI.Efficiency ab = 0.92 "Absorptivity" annotation(Dialog(group="Technical data"));

	parameter Boolean const_alpha = true "If true then constant convective heat transfer coefficient";
	parameter SI.CoefficientOfHeatTransfer alpha=1 if const_alpha "Convective heat transfer coefficient";
	// FIXME what happens if const_alpha==0?

	Medium.BaseProperties medium;

	SI.SpecificEnthalpy h_in "Specific enthalpy at inlet";
	SI.SpecificEnthalpy h_out(start=h_0) "Specific enthalpy at outlet";

	SI.Temperature T_in=Medium.temperature(state_in) "Temperature at inlet";
	SI.Temperature T_out=Medium.temperature(state_out) "Temperature at outlet";

	SI.HeatFlowRate Q_loss "Total losses";
	SI.HeatFlowRate Q_rad "Radiative losses";
	SI.HeatFlowRate Q_con "Convective losses";
	SI.HeatFlowRate Q_rcv "Heat flow captured by curtain";

	Modelica.Blocks.Interfaces.RealInput Tamb annotation (Placement(
		transformation(
		extent={{-12,-12},{12,12}},
		rotation=-90,
		origin={0,84}), iconTransformation(
		extent={{-6,-6},{6,6}},
		rotation=-90,
		origin={0,78})));

	Modelica.Blocks.Interfaces.BooleanInput on annotation (Placement(
		transformation(extent={{-38,-94},{2,-54}}), iconTransformation(extent={{
		-24,-98},{-12,-86}})));

	Real eff;

	// Solid particle geometry
	parameter SI.Length d_p = 280e-6 "particle diameter [m]" annotation(Dialog(group="Technical data"));
	parameter SI.Area a = pi*(d_p/2)^2 "cross sectional particle area [m2]";

	// Medium properties
	parameter SI.Density rho = 3300 "particle density [kg/m3]";
	parameter Real phi_max = 0.6 "packing limit";
	parameter SI.SpecificHeatCapacity cp_s = 1200 "particle specific heat capacity [J/kg-K]";

	// Environmental variables
	parameter SI.Temperature T_ref = from_degC(25) "Reference temperature [K]";
	parameter SI.CoefficientOfHeatTransfer h_c = 100 "Convective heat transfer coefficient [W/m^2-K]";
	
	//Wall properties
	parameter SI.Efficiency eps_w = 0.8 "receiver wall emissivity";
	parameter SI.ThermalConductivity k_w = 0.2 "b wall thermal conductivity [W/m-K]";
	parameter SI.Length th_w = 0.05 "b wall thickness [m]";

	// Receiver geometry
	parameter Real AR = 1 "Receiver aspect ratio";
	
	//Discretisation
	parameter Integer Nx = 11 "Number of height positions";
	parameter SI.Length dx = H_rcv/(Nx-1) "Vertical step size [m]";
		
	//Curtain geometry*
	Real phi[Nx-1] "Packing factor";
	SI.Length x[Nx] "Vertical position across the receiver height";
	
	//Flow properties and variables
	parameter SI.Velocity vp_in = 0.25 "Particle velocity at the inlet [m/s]";
	SI.Velocity vp[Nx-1] "Particle velocity [m/s]";
	SI.Length th_c[Nx-1] "Receiver depth";
	//Curtain radiation properties
	parameter SI.Efficiency eps_c = 0.9973435 "particle emissivity across the receiver height";
	parameter SI.Efficiency tau_c = 5.75335e-8 "particle tramittance across the receiver height";
	parameter SI.Efficiency abs_c = 0.9976785 "particle absorptance across the receiver height";

	SI.Temperature T_s[Nx-2] "Particle Temperature";
	SI.SpecificEnthalpy h_s[Nx-2] "Particle enthalpy";
	SI.Temperature T_w[Nx] "Wall Temperature accross the receiver height";
	SI.Temperature T_w_0 "Wall temperature at x=0";	
	
	SI.HeatFlux gc_f[Nx-1] "g curtain at the front";
	SI.HeatFlux gc_b[Nx-1] "g curtain at the back";
	SI.HeatFlux g_w[Nx-1] "g at the wall";
	SI.HeatFlux jc_f[Nx-1] "j curtain at the front";
	SI.HeatFlux jc_b[Nx-1] "j curtain at the back"; //

	SI.HeatFlux q_solar "Uniform solar flux [W/m2]";

protected
	parameter SI.Temperature T_0=from_degC(290) "Start value of temperature";
	parameter Medium.ThermodynamicState state_0=Medium.setState_pTX(1e5,T_0);
	parameter SI.SpecificEnthalpy h_0=Medium.specificEnthalpy(state_0);

	Medium.ThermodynamicState state_in=Medium.setState_phX(fluid_a.p,h_in);
	Medium.ThermodynamicState state_out=Medium.setState_phX(fluid_b.p,h_out);
equation
	medium.h=(h_in+h_out)/2;
	h_in=inStream(fluid_a.h_outflow);
	fluid_b.h_outflow=max(h_0,h_out);
	fluid_a.h_outflow=0;

	q_solar = heat.Q_flow/A;

	th_c_in = fluid_a.m_flow/(phi_max*rho*vp_in*W_rcv);
	for i in 1:(Nx-1) loop th_c[i] = th_c_in + 0.0087*x[i]; end for;

	x[1] = dx/2;
	for i in 2:(Nx-2) loop x[i] = x[i-1]+dx; end for;
	x[end-1] = H_rcv - dx/2;
	x[end] = H_rcv;

	T_w_0=T_w[1];
	T_w[end]=T_w[end-1];

	for i in 1:(Nx-1) loop
		if i==1 then
		phi[i] = fluid_a.m_flow/(W_rcv*rho*th_c[i]*vp[i]);
		vp[i] = (vp_in + (vp_in^2 + 4*g_n*(x[i] - 0))^0.5)/2;
		else
		phi[i] = phi[i-1]*rho*vp[i-1]*W_rcv*th_c[i-1]/(W_rcv*rho*th_c[i]*vp[i]);
		vp[i] = (vp[i-1] + (vp[i-1]^2 + 4*g_n*(x[i] - x[i-1]))^0.5)/2;
		end if;	

	if i==1 then
		0 = -(fluid_a.m_flow*h_s[i] - fluid_a.m_flow*h_in)/dx + W_rcv*(gc_f[i] - jc_f[i] + gc_b[i] - jc_b[i] - h_c*(T_s[i]-Tamb));
		0 = -k_w*((T_w[i+1]-T_w[i])/dx - (T_w[i]-T_w_0)/(dx/2))/dx - (g_w[i]-(eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]))/th_w + ((1/h_c+th_w/k_w)^(-1)*(T_w[i] - T_ref))/th_w;
	elseif i==Nx-1 then
		//Energy Balance
		0 = -(fluid_a.m_flow*h_out - fluid_a.m_flow*h_s[i-1])/dx + W_rcv*(gc_f[i] - jc_f[i] + gc_b[i] - jc_b[i] - h_c*(T_s[i]-Tamb));
		//Back Wall Energy Balance
		0 = -k_w*((T_w[i+1]-T_w[i])/(dx/2) - (T_w[i]-T_w[i-1])/dx)/dx - (g_w[i]-(eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]))/th_w + ((1/h_c+th_w/k_w)^(-1)*(T_w[i] - T_ref))/th_w;
	else
		0 = -(fluid_a.m_flow*h_s[i] - fluid_a.m_flow*h_s[i-1])/dx + W_rcv*(gc_f[i] - jc_f[i] + gc_b[i] - jc_b[i] - h_c*(T_s[i]-Tamb));
		0 = -k_w*((T_w[i+1]-T_w[i])/dx - (T_w[i]-T_w[i-1])/dx)/dx - (g_w[i]-(eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]))/th_w + ((1/h_c+th_w/k_w)^(-1)*(T_w[i] - T_ref))/th_w;
	end if;
	
	if i==(Nx-1) then	
		gc_f[i] = q_solar;
		jc_f[i] = (1-tau_c)*(eps_c*sigma*T_out^4 + (1-abs_c)*q_solar) + tau_c*gc_b[i];
		gc_b[i] = (eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]);
		jc_b[i] = (1-tau_c)*(eps_c*sigma*T_out^4 + (1-eps_c)*gc_b[i]) + tau_c*q_solar;
		g_w[i] = jc_b[i];	
	else
		//Radiation Heat Transfer
		T_s[i] = h_s[i]/cp_s + T_ref;
		gc_f[i] = q_solar;
		jc_f[i] = (1-tau_c)*(eps_c*sigma*T_s[i]^4 + (1-abs_c)*q_solar) + tau_c*gc_b[i];
		gc_b[i] = (eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]);
		jc_b[i] = (1-tau_c)*(eps_c*sigma*T_s[i]^4 + (1-eps_c)*gc_b[i]) + tau_c*q_solar;
		g_w[i] = jc_b[i];	
	end if;
	assert(T_w[i] > 0, "The wall tempeature cannot be negative in Kelvin!");
	end for;

	heat.T=medium.T;
	fluid_b.m_flow=-fluid_a.m_flow;
	fluid_a.p=medium.p;
	fluid_b.p=medium.p;

	Q_rad=A*sigma*em*(medium.T^4-Tamb^4);
	Q_con=A*alpha*(medium.T-Tamb);

	if on then
		Q_loss=-Q_rad-Q_con;
	else
		Q_loss=0;
	end if;

	Q_rcv=fluid_a.m_flow*(h_out-h_in);
	eff=max(Q_rcv, 0)/max(1,heat.Q_flow);


	annotation (Documentation(info="<html>
</html>", revisions="<html>
<ul>
<li>A. Shirazi and A. Fontalvo:<br>Released first version. </li>
</ul>
</html>"));
end ParticleReceiver_1D;
