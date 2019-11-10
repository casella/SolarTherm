within SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign;

  model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
  extends SolarTherm.Media.CO2.PropCO2;
    import SI = Modelica.SIunits;
    import CV = Modelica.SIunits.Conversions;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    parameter Real eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
    parameter SI.Diameter diam_rotor(fixed = false) "on-design diameter of the rotor";
    parameter SI.AngularVelocity N_design = 40000 * 0.104 "Design rotationnal speed in rad/s";
    parameter Real psi_des(fixed = false) "on-design adimensionned head";
    parameter Real PR = 3 "pressure ratio chosen";
    parameter SI.Power P_nom_des "nominal power at design point";
    parameter SI.ThermodynamicTemperature T_in_des = CV.from_degC(45) "chosen inlet temperature of the compressor at design point";
    parameter SI.AbsolutePressure p_high_des = CV.from_bar(250);
    parameter Real phi_opt = 0.0297035 "optimal adimensionned mass flow";
    parameter SI.AngularVelocity tipSpeed(fixed = false, start = 500) "tip Speed of the rotor";
    SI.AbsolutePressure p_out(start = p_high_des) "outlet pressure";
    MedPB.ThermodynamicState state_a "thermodynamic state at the entrance";
    MedPB.ThermodynamicState state_isen "thermodynamic state at the end of the isentropic compression";
    MedPB.ThermodynamicState state_b "thermodynamic state at the end of the real compresssion";
    SI.Power W_comp "power used for compression";
    SI.SpecificEntropy s_entrance "entropy at the entrance of the compressor";
    Real phi "adimensionned mass flow rate";
    Real psi "adimensionned head";
    SI.Density d_entrance(start = 267) "density at the inlet";
    SI.Efficiency eta_comp(start = eta_design) "isentropic efficiency of the compressor";
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {42, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
    protected
    parameter MedPB.ThermodynamicState state_in_des(p.fixed = false, h.fixed = false, h.start = 550000),state_isen_des(p.fixed = false, h.fixed = false),state_out_des(p.fixed = false, h.fixed = false);
  
    
    parameter SI.Power W_comp_des(fixed = false);
    parameter SI.MassFlowRate m_des(fixed = false, start = 3) "design mass flow rate in kg/s";
  initial equation
    2 * m_des = phi_opt * MedPB.density(state_in_des) * N_design * diam_rotor ^ 3;
    psi_des = (state_isen_des.h - state_in_des.h) / tipSpeed ^ 2;
    W_comp_des = m_des * (state_out_des.h - state_in_des.h);
  
    state_isen_des = MedPB.setState_psX(p_high_des, MedPB.specificEntropy(state_in_des));
    tipSpeed = diam_rotor * N_design / 2;
    state_out_des = MedPB.setState_phX(p_high_des,state_in_des.h + (state_isen_des.h - state_in_des.h) / eta_design);
  equation
    state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
    s_entrance = MedPB.specificEntropy(state_a);
    state_isen = MedPB.setState_psX(p_out, s_entrance);
    state_b = MedPB.setState_phX(p_out, state_a.h + (state_isen.h - state_a.h) / eta_comp);
    port_b.p = state_b.p;
    port_b.h_outflow = state_b.h;
    W_comp = port_a.m_flow * (state_b.h - state_a.h);
    port_a.h_outflow = 0;
    d_entrance = MedPB.density(state_a);
    port_a.m_flow + port_b.m_flow = 0;
    phi = port_a.m_flow / (d_entrance * tipSpeed * diam_rotor ^ 2);
    psi = (state_isen.h - state_a.h) / tipSpeed ^ 2;
    psi = (0.04049 + 54.7 * phi - 2505 * phi ^ 2 + 53224 * phi ^ 3 - 498626 * phi ^ 4) * psi_des / 0.46181921979961293;
    eta_comp = eta_design / 0.677837 * ((-0.7069) + 168.6 * phi - 8089 * phi ^ 2 + 182725 * phi ^ 3 - 1.638 * 10 ^ 6 * phi ^ 4);
    annotation(
      Diagram(graphics = {Text(origin = {-20, 18}, extent = {{-28, 16}, {42, -46}}, textString = "COMPRESSOR"), Polygon(origin = {-12, 10}, points = {{-42, 40}, {-42, -44}, {42, -70}, {42, 70}, {-42, 40}, {-42, 40}})}, coordinateSystem(initialScale = 0.1)),
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Polygon(origin = {-26, -2}, points = {{-40, 42}, {-42, -48}, {42, -78}, {42, 78}, {-40, 42}}), Text(origin = {-16, 11}, extent = {{-48, -31}, {24, 15}}, textString = "COMPRESSOR")}));
    annotation(
      Documentation(info = "<html>
  		<p>This compressor's model is based on the phD thesis of J. Dyreby.&nbsp;</p>
  <p>The performance maps comes from the Sandia National Laboratory first compressor. It should be updated. The performance maps is compressed in three correlations, expressing the adimensionned head and the efficiency as functions of the adimensionned mass flow.&nbsp;</p>
  <p>The same correlations are used; only the maximal values are changed.</p>
  <p>J. J. Dyreby, &laquo; Modeling the supercritical carbon dioxide Brayton cycle with recompression &raquo;, The University of Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
  		
  		</html>"));
  end CompressorOnShaft;