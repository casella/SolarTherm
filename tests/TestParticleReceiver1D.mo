model TestParticleReceiver1D
	import SolarTherm.{Models,Media};
	import Modelica.SIunits.Conversions.from_degC;
	import SI = Modelica.SIunits;
	import nSI = Modelica.SIunits.Conversions.NonSIunits;
	import CN = Modelica.Constants;
	import CV = Modelica.SIunits.Conversions;
	import FI = SolarTherm.Models.Analysis.Finances;
	import Modelica.Math;
	import Modelica.Blocks;
	replaceable package Medium = SolarTherm.Media.SolidParticles.CarboHSP_ph "Medium props for Carbo HSP 40/70";
	
	parameter SI.HeatFlowRate Q_in = 1200 * 788.8 * 593.7 "Incident heat flux on receiver surface, in W";
	parameter SI.MassFlowRate m_in = 1827 "Receiver inlet mass flow rate, in kg/s";
	
	Modelica.Fluid.Sources.FixedBoundary source(
		redeclare package Medium = Medium, 
		T = Modelica.SIunits.Conversions.from_degC(580.3), 
		nPorts = 1, 
		p = 1e5, 
		use_p = true) 
		annotation(
		Placement(visible = true, transformation(origin = {58, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));

	Modelica.Fluid.Sources.FixedBoundary sink(
		redeclare package Medium = Medium, 
		nPorts = 1, 
		p = 1e5) 
		annotation(
		Placement(visible = true, transformation(extent = {{32, 22}, {12, 42}}, rotation = 0)));

	SolarTherm.Models.CSP.CRS.Receivers.ParticleReceiver_1D receiver(
		redeclare package Medium = Medium) 
		annotation(
		Placement(visible = true, transformation(extent = {{-40, 4}, {-4, 40}}, rotation = 0)));

	Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow Heat(
		Q_flow = Q_in) 
		annotation(
		Placement(visible = true, transformation(origin = {-80, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Sources.BooleanExpression Operation(
		y = true) 
		annotation(
		Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Sources.RealExpression Tamb(
		y = 25) 
		annotation(
		Placement(visible = true, transformation(origin = {-46, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	SolarTherm.Models.Fluid.Pumps.ParticleLift particleLift(
		m_flow_fixed = 1827, 
		use_input = false) 
		annotation(
		Placement(visible = true, transformation(origin = {10, -14}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));

equation
	connect(particleLift.port_a, source.ports[1]) 
		annotation(
		Line(points = {{20, -14}, {48, -14}}, color = {0, 127, 255}));
	connect(particleLift.port_b, receiver.fluid_a) 
		annotation(
		Line(points = {{0, -14}, {-18, -14}, {-18, 6}}, color = {0, 127, 255}));
	connect(sink.ports[1], receiver.fluid_b) 
		annotation(
		Line(points = {{12, 32}, {-16, 32}, {-16, 31}}, color = {0, 127, 255}));
	connect(receiver.heat, Heat.port) 
		annotation(
		Line(points = {{-40, 27}, {-57, 27}, {-57, 44}, {-70, 44}}, color = {191, 0, 0}));
	connect(Operation.y, receiver.on) 
		annotation(
		Line(points = {{-69, -4}, {-39.5, -4}, {-39.5, 5}, {-25, 5}}, color = {255, 0, 255}));
	connect(Tamb.y, receiver.Tamb) 
		annotation(
		Line(points = {{-34, 66}, {-34, 51}, {-22, 51}, {-22, 36}}, color = {0, 0, 127}));
	annotation(
		uses(Modelica(version = "3.2.2"), SolarTherm(version = "0.2")));
end TestParticleReceiver1D;
