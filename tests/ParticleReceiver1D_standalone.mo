
function function_1 "function used in curtain property calculation"
  input Real eps_s "emissivity of the particle";
  input Real Cta "product C[i]*th_c[i]*a";
  output Real result;
algorithm
  result := 1-(1/((eps_s*Cta)^2))+((1+eps_s*Cta)/((eps_s*Cta)^2))*exp(-eps_s*Cta);
end function_1;

model ParticleReceiver1D_standalone "Falling particle flow and energy model"
	/* FIXME this model uses fixed c_p in the finite different equations, but
	   nonlinear h(T) for the inlet/outlet, leading to inconsistencies. But it
	   solves! */

	/* FIXME model is lacking connectors etc so is not yet ready for integration
	   into larger models. */

	import SI = Modelica.SIunits;
	import CONST = Modelica.Constants;
	import Modelica.SIunits.Conversions.*;
	import SolarTherm.Media;

	constant Boolean fixed_geometry = false;
    constant Boolean with_wall_conduction = false;
	constant Boolean fixed_cp = false;
	constant Boolean with_isothermal_backwall = false;
	constant Boolean with_uniform_curtain_props = false;
	constant SI.SpecificHeatCapacity cp_s = 1200. "solid specific heat capacity [J/kg-K]";

	//Discretisation
	parameter Integer N = 20 "Number of vertical elements";

	// Medium
	replaceable package Medium = Media.SolidParticles.CarboHSP_ph;

	// temperature used to initialise screen
	parameter SI.Temperature T_ref = from_degC(580.3);
	parameter SI.SpecificEnthalpy h_0 = (
			if fixed_cp then 0
			else Medium.specificEnthalpy(Medium.setState_pTX(p_des,T_ref))
		);

	// Solid particle geometry
	parameter SI.Length d_p = 280e-6 "particle diameter [m]" annotation(Dialog(group="Technical data"));
	parameter SI.Area a = 0.25 * CONST.pi * d_p^2 "cross sectional particle area [m2]";

	// Medium properties
	parameter SI.Efficiency eps_s = 0.86 "solid particle emissivity";
	parameter SI.Efficiency abs_s = 0.92 "solid particle absorptivity";
	parameter SI.Density rho_s = 3300. "solid density [kg/m3]";
	parameter Real phi_s_max = 0.6 "packing limit";

	// Environmental variables
	parameter SI.Temperature T_amb = from_degC(25) "Ambient temperature [K]";
	parameter SI.CoefficientOfHeatTransfer h_conv = 100. "Convective heat transfer coefficient [W/m^2-K]";

	//Wall properties
	parameter SI.Efficiency eps_w = 0.8 "receiver wall emissivity";
	parameter SI.ThermalConductivity k_w = 0.2 "b wall thermal conductivity [W/m-K]";
	parameter SI.Length t_w = 0.05 "back wall thickness [m]";

	// Receiver geometry
	parameter Real AR = 1 "Receiver aspect ratio";
	//parameter SI.Length t_c_in = 1887.76/(0.6*1200*0.25*24.37) "Curtain thicknesss at the inlet";
	SI.Length t_c_in (start=1887.76/(0.6*1200*0.25*24.37)) "Curtain thicknesss at the inlet";

	SI.Length H_drop "Receiver drop height [m]";
	SI.Area A_ap "Receiver aperture area [m2]";
	SI.Length w_c "Aperture (curtain) width [m], w_curtain";
	SI.Length dx "Vertical step size [m]";

	// fixed for now, later it would be variable...
	parameter SI.Temperature T_in = from_degC(580.3) "Inlet temperature [K]";
	//parameter SI.Temperature T_in = from_degC(150) "Inlet temperature [K]";
	//parameter SI.Temperature T_in = T_amb  "Inlet temperature [K]";
	parameter SI.Temperature T_out = from_degC(800.) "Outlet temperature [K]";
	//parameter SI.Temperature T_in_des = from_degC(580.3) "Design inlet temperature, K";
	parameter SI.Pressure p_des = 1e5 "Design pressure, Pa";

	//Variables
	SI.MassFlowRate mdot(start=1600,nominal=1000) "Inlet mass flow rate [kg/s]";

	//Curtain geometry
	Real phi_s__[N+1] (start=fill(0.5,N+1),min=zeros(N+1),max=fill(1,N+1)) "Curtain packing factor";

	// Vertical positions
	SI.Length x__[N+2] (min=zeros(N+2),max=fill(100.,N+2)) "Vertical positions of nodes";

	// NOTE variables with a '__' appended to the name needed a [0] index in EES. Hence they are index with [i+1] in this code.

	//Flow properties and variables
	parameter SI.Velocity v_s_in = 0.25 "Inlet curtain velocity [m/s]";
	SI.Velocity v_s__[N+1] (start=fill(1.5*v_s_in,N+1),min=fill(v_s_in,N+1),max=fill(1000,N+1)) "Particles velocity [m/s]";
	SI.Length t_c__[N+2] "Receiver depth";
	//Real C[N] "FIXME something to do with curtain opacity";
	SI.Temperature T_s__[N+1] (start = fill(T_in,N+1), max=fill(3000.,N+1)) "Curtain Temperature";
	SI.SpecificEnthalpy h_s__[N+1] (start = fill(h_0,N+1), max=fill(cp_s*(2000-T_ref),N+1)) "Curtain enthalpy";
	SI.SpecificEnthalpy h_out (start = h_0) "Curtain outlet enthalpy";

	//Wall properties and variables
	SI.Temperature T_w__[N+2] (start = fill(T_amb+1.,N+2), max=fill(3000.,N+2)) "Receiver wall temperature";
	//SI.Temperature T_w_0 (start = T_in) "Wall inlet temperature";

	//Curtain radiation properties
	SI.Efficiency eps_c[N] (start=fill(0.5,N),max=fill(1.,N),min=fill(0.,N)) "Curtain emissivity";
	SI.Efficiency tau_c[N] (start=fill(0.5,N),max=fill(1.,N),min=fill(0.,N)) "Curtain tramittance";
	SI.Efficiency abs_c[N] (start=fill(0.5,N),max=fill(1.,N),min=fill(0.,N)) "Curtain absorptance";

	//Radiation heat fluxes
	SI.HeatFlux q_solar "Uniform solar flux [W/m2]";

	//Overall performance
	SI.HeatFlowRate Qdot_rec "Total heat rate absorbed by the receiver";
	SI.HeatFlowRate Qdot_inc "Total heat rate incident upon the receiver (before losses)";
	Real eta_rec(min=0, max=1) "Receiver efficiency";

//protected

	SI.HeatFlux gc_f[N] (min=zeros(N)) "Curtain radiation gain at the front";
	SI.HeatFlux jc_f[N] (min=zeros(N)) "Curtain radiation loss at the front";
	SI.HeatFlux gc_b[N] (min=zeros(N)) "Curtain radiation gain at the back";
	SI.HeatFlux jc_b[N] (min=zeros(N)) "Curtain radiation loss at the back";
	SI.HeatFlux g_w[N] (min=zeros(N)) "Wall radiation gain from the front";
	SI.HeatFlux j_w[N] (min=zeros(N)) "Wall radiosity (front)";

	SI.HeatFlux q_conv[N] "Heat flux loss by convection";
	SI.HeatFlux q_net_c[N] "Net heat gain by curtain";

	SI.MassFlowRate mdot_check "mass balance check";
	SI.HeatFlowRate Qdot_check "heat balance check";
equation
	dx * N = H_drop;
	AR * H_drop = w_c;
	A_ap = H_drop * w_c;

	//q_solar = 1200*788.8;///A_ap; //q_solar=DNI*CR
	q_solar = 1200 * 788.8;

	h_out = h_s__[N]; // just an alias

	if fixed_geometry then
		// specify the geometry and flow rate, see how the particles are heated
		H_drop = 24.37;
		mdot = 1887.76;
	else
		// specify the final temperature and heat flux.
		mdot = 1827.;
		T_s__[N+1] = T_out;
	end if;

	// Boundary conditions
	phi_s__[1] = phi_s_max;
	//t_c__[1] = t_c_in; // covered by loop below
	//phi_s__[1]*rho_s*v_s__[1]*t_c__[1]*w_c = mdot; // to determine th_c[1]
	T_s__[1] = T_in;
	// h_s__[1] = Medium.specificEnthalpy(Medium.setState_pTX(p_des,T_s__[1])); // covered in loop below
	v_s__[1] = v_s_in;
	T_w__[1] = T_w__[2];
	T_w__[N+2]=T_w__[N+1];

	x__[1] = 0;
	for i in 2:N+1 loop
		x__[i] = dx*(1./2 + (i-2));
	end for;
	x__[N+2] = H_drop;

	for i in 1:N+1 loop
		// Mass balance
		mdot = phi_s__[i]*rho_s*v_s__[i]*t_c__[i]*w_c;
	end for;

	for i in 2:N+1 loop
		// Curtain Mass balance
		//phi_s__[i]*t_c__[i]*rho_s*v_s__[i]
		//	= phi_s__[i-1]*t_c__[i-1]*rho_s*v_s__[i-1];

		// Curtain momentum balance (curtain opacity)
		mdot*(v_s__[i] - v_s__[i-1])
				= dx*w_c*phi_s__[i]*t_c__[i]*rho_s* CONST.g_n;

	end for;

	for i in 1:N+1 loop
		h_s__[i] =
		 	if fixed_cp then cp_s*(T_s__[i]-T_ref)
			else Medium.specificEnthalpy(Medium.setState_pTX(p_des,T_s__[i]));
	end for;

	for i in 1:N+2 loop
		// Curtain thickness, FIXME note magic number 'growth factor'...
		t_c__[i] = t_c_in + 0.0087*x__[i];
	end for;

	for i in 1:N loop

		if with_uniform_curtain_props then
			eps_c[i] = eps_s;
			abs_c[i] = abs_s;
			tau_c[i] = 0.4; //exp(-3*phi_s__[i+1]*t_c__[i+1]/(2*d_p));
		else
			// Curtain radiation properties
			eps_c[i]*(1-tau_c[i]) = function_1(eps_s, 6*phi_s__[i+1]/(CONST.pi*d_p^3) * t_c__[i+1] * a);
			abs_c[i]*(1-tau_c[i]) = function_1(abs_s, 6*phi_s__[i+1]/(CONST.pi*d_p^3) * t_c__[i+1] * a);
			tau_c[i] = exp(-3*phi_s__[i+1]*t_c__[i+1]/(2*d_p));
		end if;

		// Curtain energy balance
		q_net_c[i] =  gc_f[i] - jc_f[i] + gc_b[i] - jc_b[i] - h_conv*(T_s__[i+1]-T_amb);

		//q_net_c[i]*dx = phi_s__[i+1]*t_c__[i+1]*rho_s*v_s__[i+1]*h_s__[i+1] - phi_s__[i]*t_c__[i]*rho_s*v_s__[i]*h_s__[i];
		q_net_c[i]*dx*w_c = mdot*(h_s__[i+1] - h_s__[i]);

		// Curtain-wall radiation heat fluxes (W/m²)
		gc_f[i] = q_solar;
		jc_f[i] = (1-tau_c[i])*(eps_c[i] * CONST.sigma * T_s__[i+1]^4 + (1-abs_c[i])*q_solar) + tau_c[i]*gc_b[i];
		gc_b[i] = j_w[i];
		jc_b[i] = (1-tau_c[i])*(eps_c[i] * CONST.sigma * T_s__[i+1]^4 + (1-eps_c[i])*gc_b[i]) + tau_c[i]*q_solar;
		g_w[i] = jc_b[i];
		j_w[i] = eps_w * CONST.sigma * T_w__[i+1]^4 + (1-eps_w)*g_w[i];

		// Back wall energy balance
		if with_isothermal_backwall then
			// wall is at ambient temperature, absorbed heat lost as convection+radiation
			T_w__[i+1] = T_amb;
			q_conv[i] + j_w[i] = g_w[i];
		else
			q_conv[i] = (T_w__[i+1] - T_amb) / (1/h_conv + t_w/k_w);
			0 = (
					if with_wall_conduction
					then -k_w*(((T_w__[i+2]-T_w__[i+1])/(x__[i+2]-x__[i+1])) - ((T_w__[i+1]-T_w__[i])/(x__[i+1]-x__[i])))*t_w
					else 0
				) - (g_w[i]-(eps_w * CONST.sigma * T_w__[i+1]^4 + (1-eps_w)*g_w[i]))*dx
				+ (q_conv[i])*dx;
		end if;
	end for;

	mdot_check = phi_s__[1]*rho_s*v_s__[1]*w_c*t_c__[1]  -  phi_s__[N+1]*rho_s*v_s__[N+1]*w_c*t_c__[N+1];

	Qdot_inc = q_solar * A_ap;
	Qdot_rec = mdot * (h_s__[N+1] - h_s__[1]);
	eta_rec = Qdot_rec / Qdot_inc;

	Qdot_check = Qdot_rec - sum(dx*w_c*q_net_c[i] for i in 1:N);

	annotation (Documentation(info="<html>
<p>Model based on an EES-based model written by Kevin Albrecht at Sandia
National Laboratories.</p></html>", revisions="<html>
<dl>
<dt>Armando Fontalvo:<dd>Initial development as a stand-alone model.
<dt>John Pye:<dd>Expanded model with various binary switches plus variable geometry mode.
</ul>
</html>"));
end ParticleReceiver1D_standalone;
