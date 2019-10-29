within SolarTherm.Models.Fluid.HeatExchangers;
model SimpleHeatExchanger "A simple counterflow heat exchanger model based on LMTD method"
	extends SolarTherm.Interfaces.Models.HeatExchangerFluid;
	//extends SolarTherm.Interfaces.Models.HeatExchangerFluid(port_b_out.h_outflow(start=1200e3));

	import CN = Modelica.Constants;
	import SI = Modelica.SIunits;
	import Modelica.Math;

	parameter SI.Area A = 1 "Heat transfer surface area";
	parameter SI.CoefficientOfHeatTransfer U = 1 "Heat tranfer coefficient";
	parameter SI.Temperature dT_approach = 1 "Approach temperature";

	Medium_A.BaseProperties medium_a_in;
	Medium_A.BaseProperties medium_a_out;
	Medium_B.BaseProperties medium_b_in;
	Medium_B.BaseProperties medium_b_out;

	SI.Temperature T_a_in(start=800+273.15) "Medium A inlet temperature";
	SI.Temperature T_a_out(start=580.3+273.15) "Medium A outlet temperature";
	SI.Temperature T_b_in(start=565.3+273.15) "Medium B inlet temperature";
	SI.Temperature T_b_out(start=715+273.15) "Medium B outlet temperature";

	SI.HeatFlowRate Q_flow "Heat flow from hot to cold side";
	SI.Temperature LMTD(start=40.36) "Logarithmic mean temperature difference";

equation
	port_a_in.m_flow + port_a_out.m_flow = 0;
	port_b_in.m_flow + port_b_out.m_flow = 0;

	port_a_out.Xi_outflow = inStream(port_a_in.Xi_outflow);
	port_a_in.Xi_outflow = inStream(port_a_out.Xi_outflow);
	port_b_out.Xi_outflow = inStream(port_b_in.Xi_outflow);
	port_b_in.Xi_outflow = inStream(port_b_out.Xi_outflow);

	port_a_out.C_outflow = inStream(port_a_in.C_outflow);
	port_a_in.C_outflow = inStream(port_a_out.C_outflow);
	port_b_out.C_outflow = inStream(port_b_in.C_outflow);
	port_b_in.C_outflow = inStream(port_b_out.C_outflow);

	medium_a_in.p = port_a_in.p;
	medium_a_out.p = port_a_out.p;
	medium_b_in.p = port_b_in.p;
	medium_b_out.p = port_b_out.p;

	medium_a_in.h = inStream(port_a_in.h_outflow);
	medium_a_out.h = port_a_out.h_outflow;
	medium_b_in.h = inStream(port_b_in.h_outflow);
	medium_b_out.h = port_b_out.h_outflow;

	medium_a_in.T = T_a_in;
	medium_a_out.T = T_a_out;
	medium_b_in.T = T_b_in;
	medium_b_out.T = T_b_out;

	dT_approach = T_a_out - T_b_in;
	Q_flow = port_a_in.m_flow*(inStream(port_a_in.h_outflow) - port_a_out.h_outflow);
	Q_flow = U*A*LMTD;
	LMTD = ((T_a_in-T_b_out)-(T_a_out-T_b_in))/(Math.log((T_a_in-T_b_out)/(T_a_out-T_b_in)));
	Q_flow = port_b_in.m_flow*(port_b_out.h_outflow - inStream(port_b_in.h_outflow));

	port_a_out.p = port_a_in.p;
	port_b_out.p = port_b_in.p;

	// Shouldn't have reverse flows
	port_a_in.h_outflow = 0.0;
	port_b_in.h_outflow = 0.0;

	annotation (Documentation(info="<html>
</html>", revisions="<html>
<ul>
<li>A. Shirazi:<br>Released first version. </li>
</ul>
</html>"));

end SimpleHeatExchanger;
