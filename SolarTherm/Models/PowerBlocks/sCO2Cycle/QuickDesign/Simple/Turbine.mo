within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple;

 model Turbine "OD model of a turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter SI.Efficiency eta_design = 0.9 "isentropic efficiency of the turbine";
      parameter SI.Efficiency PR = 3 "Pressure ratio";
      parameter SI.AngularVelocity N_shaft = 3358;
      parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "TIT at design point";
      //parameters for off and on-design
      parameter Integer N_points = 10;
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false, each start = 2900) "inlet entropy to the turbine";
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false);
      parameter SI.Power[N_points + 1] W_turb(each fixed = false);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false);
    initial equation
//Part I. On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p / PR, state_in[1].h + eta_design * (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[1].p / PR, s_in[1])) - state_in[1].h));
      W_turb[1] = m_des * (state_out[1].h - state_in[1].h);
      for i in 2:N_points + 1 loop
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_out[i] = MedPB.setState_phX(state_in[i].p / PR, state_in[i].h + eta_design * (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[i].p / PR, s_in[i])) - state_in[i].h));
        W_turb[i] = m_des * (state_out[i].h - state_in[i].h);
      end for;
      annotation(
        Documentation(info = "<html>
  		<p>This turbine's model is based on the phD thesis of J. Dyreby.&nbsp;</p>
  <p>The isentropic efficiency is calculated as a function of the tip speed ration between the tip speed of the rotor and the spouting velocity. It is said to be functionnal for any size.</p>
  <p>The outlet pressure goes beyond the critical pressure for a mass flow too small. The cycle calculation should therefore not be performed below this pressure.</p>
  <p>J. J. Dyreby, &laquo; Modeling the supercritical carbon dioxide Brayton cycle with recompression &raquo;, The University of Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
  		</html>"));
      annotation(
        Diagram(graphics = {Text(origin = {-36, -28}, extent = {{18, 80}, {78, 16}}, textString = "TURBINE"), Polygon(origin = {15, 20}, points = {{-35, 44}, {-35, -52}, {35, -68}, {35, 68}, {-35, 44}, {35, 68}, {-35, 44}})}, coordinateSystem(initialScale = 0.1)),
  Icon(graphics = {Text(origin = {-10, 26}, extent = {{-10, 12}, {52, -34}}, textString = "TURBINE"), Ellipse(extent = {{56, 58}, {56, 58}}, endAngle = 360), Polygon(origin = {11, 17}, points = {{-37, 49}, {-37, -51}, {37, -71}, {37, 71}, {-37, 49}})}, coordinateSystem(initialScale = 0.1)));
    end Turbine;