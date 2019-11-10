within SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign;

model Turbine "OD model of a turbine"
    extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    import SI = Modelica.SIunits;
    parameter SI.Efficiency eta_design = 0.9 "isentropic efficiency of the turbine";
    parameter SI.Efficiency PR = 3 "Pressure ratio";
  
    parameter Modelica.SIunits.Area A_nozzle(fixed = false);
    parameter SI.AngularVelocity N_shaft = 3358;
    parameter SI.Diameter diam_turb(fixed = false);
    parameter SI.Velocity tipSpeed_des(fixed = false);
  
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {-32, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-38, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {60, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SI.Efficiency eta_turb "efficiency of the turbine";
      SI.Density d_outlet;
      SI.Velocity C_spouting(start = C_spouting_des);
      MedPB.ThermodynamicState state_a (h.start=state_in_des.h) "thermodynamic state at the entrance";
      MedPB.ThermodynamicState state_isen "thermodynamic state at the end of the isentropic decompression";
      MedPB.ThermodynamicState state_b "thermodynamic state at the end of the real decompresssion";
      SI.Power W_turb "Outlet power";
      SI.AbsolutePressure p_out(start = state_out_des.p);
      SI.SpecificEntropy s_entrance " entropy at the entrance of the turbine";
      
    protected
    parameter MedPB.ThermodynamicState state_in_des(p.fixed = false, h.fixed = false, h.start = 1.2 * 10 ^ 6), state_isen_des(p.fixed = false, h.fixed = false), state_out_des(p.fixed = false, h.fixed = false, h.start = 900000);
  
    parameter SI.Velocity C_spouting_des(fixed = false, start = 500);
    parameter SI.MassFlowRate m_des(fixed = false);
    parameter SI.Power W_turb_des(fixed = false);
  initial equation
    state_isen_des = MedPB.setState_psX(state_in_des.p / PR, MedPB.specificEntropy(state_in_des));
    state_out_des = MedPB.setState_phX(state_in_des.p / PR,state_in_des.h + (state_isen_des.h - state_in_des.h) * eta_design);
    C_spouting_des ^ 2 = 2 * (state_in_des.h - state_isen_des.h);
    m_des = C_spouting_des * A_nozzle * MedPB.density(state_out_des);
    W_turb_des = m_des * (state_out_des.h - state_in_des.h);
    tipSpeed_des = N_shaft * diam_turb / 2;
    tipSpeed_des / C_spouting_des = 0.707;
  equation
  state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
  s_entrance = MedPB.specificEntropy(state_a);
  state_isen = MedPB.setState_psX(p_out, s_entrance);
  state_b = MedPB.setState_phX(p_out, state_a.h + (state_isen.h - state_a.h) * eta_turb);
  port_b.p = state_b.p;
  
  port_a.h_outflow = inStream(port_b.h_outflow);
  port_b.h_outflow = state_b.h;
  W_turb = port_a.m_flow * (state_b.h - state_a.h);
  port_a.m_flow + port_b.m_flow = 0;
  d_outlet = MedPB.density(state_b);
  port_a.m_flow = C_spouting * A_nozzle * d_outlet;
  C_spouting ^ 2 = 2 * (state_a.h - state_isen.h);
  eta_turb = eta_design * 2 * (tipSpeed_des / C_spouting) * sqrt(1 - (tipSpeed_des / C_spouting) ^ 2);
  
    annotation(
      Diagram(graphics = {Text(origin = {-36, -28}, extent = {{18, 80}, {78, 16}}, textString = "TURBINE"), Polygon(origin = {15, 20}, points = {{-35, 44}, {-35, -52}, {35, -68}, {35, 68}, {-35, 44}, {35, 68}, {-35, 44}})}, coordinateSystem(initialScale = 0.1)),
  Icon(graphics = {Text(origin = {-10, 26}, extent = {{-10, 12}, {52, -34}}, textString = "TURBINE"), Ellipse(extent = {{56, 58}, {56, 58}}, endAngle = 360), Polygon(origin = {11, 17}, points = {{-37, 49}, {-37, -51}, {37, -71}, {37, 71}, {-37, 49}})}, coordinateSystem(initialScale = 0.1)));
    annotation(
      Documentation(info = "<html>
  		<p>This turbine's model is based on the phD thesis of J. Dyreby.&nbsp;</p>
  <p>The isentropic efficiency is calculated as a function of the tip speed ration between the tip speed of the rotor and the spouting velocity. It is said to be functionnal for any size.</p>
  <p>The outlet pressure goes beyond the critical pressure for a mass flow too small. The cycle calculation should therefore not be performed below this pressure.</p>
  <p>J. J. Dyreby, &laquo; Modeling the supercritical carbon dioxide Brayton cycle with recompression &raquo;, The University of Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
  		</html>"));
  end Turbine;