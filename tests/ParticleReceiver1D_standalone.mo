model ParticleReceiver1D_standalone "Falling particle flow and energy model"
	import SI = Modelica.SIunits;
	import Modelica.Constants.*;
	import Modelica.SIunits.Conversions.*;
	import SolarTherm.Media;
	
	/* FIXME this model uses fixed c_p in the finite different equations, but
	   nonlinear h(T) for the inlet/outlet, leading to inconsistencies. But it
	   solves! */

	/* FIXME model is lacking connectors etc so is not yet ready for integration
	   into larger models. */

	replaceable package Medium = Media.SolidParticles.CarboHSP_ph;
	// Use of medium
	Medium.ThermodynamicState state_a = Medium.setState_pTX(1e5,from_degC(580.3));
	SI.SpecificEnthalpy h0 = Medium.specificEnthalpy(state_a);
	
	// Solid particle geometry
	parameter SI.Length d_p = 280e-6 "particle diameter [m]" annotation(Dialog(group="Technical data"));
	parameter SI.Area a = pi*(d_p/2)^2 "cross sectional particle area [m2]";

	// Medium properties
	parameter SI.Efficiency eps_s = 0.86 "solid particle emissivity";
	parameter SI.Efficiency abs_s = 0.92 "solid particle absorptivity";
	parameter SI.Density rho_s = 3300 "solid density [kg/m3]";
	parameter Real phi_s_max = 0.6 "packing limit";
	parameter SI.SpecificHeatCapacity cp_s = 1200 "solid specific heat capacity [J/kg-K]";

	// Environmental variables
	parameter SI.Temperature T_0 = 298.15 "Reference temperature [K]";
	parameter SI.CoefficientOfHeatTransfer h_c = 100 "Convective heat transfer coefficient [W/m^2-K]";
	
	//Wall properties
	parameter SI.Efficiency eps_w = 0.8 "receiver wall emissivity";
	parameter SI.ThermalConductivity k_w = 0.2 "b wall thermal conductivity [W/m-K]";
	parameter SI.Length th_w = 0.05 "b wall thickness [m]";

	// Receiver geometry
	parameter Real AR = 1 "Receiver aspect ratio";
	parameter SI.Length w_curtain = 24.37 "Aperture width [m]";
	parameter SI.Length H_drop = w_curtain/AR "Receiver drop height [m]";
	parameter SI.Area A_ap = H_drop*w_curtain "Receiver aperture area [m2]";
	parameter SI.Length th_c_in = 1887.76/(0.6*1200*0.25*24.37) "Curtain thicknesss at the inlet";
	
	//Discretisation
	parameter Integer N = 100 "Number of height positions";
	parameter SI.Length dx = H_drop/(N) "Vertical step size [m]";
	
	//Variables
	SI.MassFlowRate m_flow "Inlet mass flow rate [kg/s]";
	SI.Temperature T_a "Inlet temperature [K]";
			
	//Curtain geometry*
	Real phi_s[N] "Packing factor";
	
	SI.Length x[N+1] "Vertical position across the receiver height";
	
	//Flow properties and variables
	SI.Velocity v_s_in "Inlet solid particles velocity [m/s]";
	SI.Velocity v_s[N] "Particles velocity [m/s]";
	SI.Length th_c[N] "Receiver depth";
	SI.Temperature T_s[N] (start = fill(from_degC(580.3),N)) "Particle Temperature";
	SI.SpecificEnthalpy h_s[N] (start = fill(cp_s*(580.3 - 298.15),N)) "Particle enthalpy";
	SI.SpecificEnthalpy h_in (start = cp_s*(580.3 - 298.15)) "Inlet particle enthalpy";
	SI.SpecificEnthalpy h_out (start = cp_s*(580.3 - 298.15)) "Outlet particle enthalpy";
	
	//Wall properties and variables
	SI.Temperature T_w[N+1] (start = fill(from_degC(580.3),N+1)) "Wall Temperature accross the receiver height";
	SI.Temperature T_w_0 (start = from_degC(580.3)) "Wall temperature at x=0";
	
	//Curtain radiation properties
	SI.Efficiency eps_c[N] "solid particle emissivity across the receiver height";
	SI.Efficiency tau_c[N] "solid particle tramittance across the receiver height";
	SI.Efficiency abs_c[N] "solid particle absorptance across the receiver height";
	
	//Radiation heat fluxes
	SI.HeatFlux q_solar "Uniform solar flux [W/m2]";
	
protected
	SI.HeatFlux gc_f[N] "g curtain at the front";
	SI.HeatFlux gc_b[N] "g curtain at the back";
	SI.HeatFlux g_w[N] "g at the wall";
	SI.HeatFlux jc_f[N] "j curtain at the front";
	SI.HeatFlux jc_b[N] "j curtain at the back";	 
equation
		
	q_solar = 1200*788.8;///A_ap; //q_solar=DNI*CR
	v_s_in = 0.25;
	h_in = cp_s*(T_a-T_0);	
	h_out = h_s[N];
	
	m_flow = 1887.76;
	T_a = from_degC(580.3);
	
	x[1] = dx/2;
	for i in 2:(N-1) loop x[i] = x[i-1]+dx; end for;
	x[end-1] = H_drop - dx/2;
	x[end] = H_drop;	
	
	T_w_0=T_w[1];
	T_w[end]=T_w[end-1];
	
	for i in 1:(N) loop
      	0 = - phi_s[i] + (if i==1 then (phi_s_max*th_c_in*v_s_in)/(th_c[i]*v_s[i]) else (phi_s[i-1]*th_c[i-1]*v_s[i-1])/(th_c[i]*v_s[i]));
		0 = - v_s[i] + (if i==1 then ((v_s_in + (v_s_in^2 + 4*g_n*(x[i] - 0))^0.5)/2) else (v_s[i-1] + (v_s[i-1]^2 + 4*g_n*(x[i] - x[i-1]))^0.5)/2);
	
		th_c[i] = th_c_in + 0.0087*x[i];
		
		eps_c[i]*(1-tau_c[i]) = 1-(1/((eps_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)^2))+((1+eps_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)/((eps_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)^2))*exp(-eps_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a);
		
		abs_c[i]*(1-tau_c[i]) = 1-(1/((abs_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)^2))+((1+abs_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)/((abs_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a)^2))*exp(-abs_s*(6*phi_s[i]/(pi*d_p^3))*(th_c[i])*a);
		
		tau_c[i] = exp(-3*phi_s[i]*th_c[i]/(2*d_p));
		
		//Energy Balance
		0 = -(phi_s[i]*th_c[i]*rho_s*v_s[i]*h_s[i] - (if i==1 then phi_s_max*th_c_in*rho_s*v_s_in*h_in else phi_s[i-1]*th_c[i-1]*rho_s*v_s[i-1]*h_s[i-1]))/dx + gc_f[i] - jc_f[i] + gc_b[i] - jc_b[i] - h_c*(T_s[i]-T_0);
				
		h_s[i] = cp_s*(T_s[i]-T_0);
		
		//Radiation Heat Transfer
		gc_f[i] = q_solar;
		jc_f[i] = (1-tau_c[i])*(eps_c[i]*sigma*T_s[i]^4 + (1-abs_c[i])*q_solar) + tau_c[i]*gc_b[i];
		gc_b[i] = (eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]);
		jc_b[i] = (1-tau_c[i])*(eps_c[i]*sigma*T_s[i]^4 + (1-eps_c[i])*gc_b[i]) + tau_c[i]*q_solar;
		g_w[i] = jc_b[i];
		
		//Back Wall Energy Balance
		0 = -k_w*((if i==N then 0 else (T_w[i+1]-T_w[i])/(x[i+1]-x[i])) - (if i==1 then 0 else (T_w[i]-T_w[i-1])/(x[i]-x[i-1])))/dx - (g_w[i]-(eps_w*sigma*T_w[i]^4 + (1-eps_w)*g_w[i]))/th_w + ((1/h_c+th_w/k_w)^(-1)*(T_w[i] - T_0))/th_w;
	end for;

	annotation (Documentation(info="<html>
<p>Model based on an EES-based model written by Kevin Albrecht at Sandia
National Laboratories.</p></html>", revisions="<html>
<ul>
<li>Armando Fontalve: <br>Initial development as a stand-alone model. </li>
</ul>
</html>"));
end ParticleReceiver1D_standalone;

