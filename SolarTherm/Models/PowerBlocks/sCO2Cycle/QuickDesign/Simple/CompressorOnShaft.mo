within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple;

 model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      parameter SI.Efficiency eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
      parameter SI.AngularVelocity N_design = 40000 * 0.104 "Design rotationnal speed in rad/s";
      parameter Real PR = 3 "pressure ratio chosen";
      //On-design calculations
      parameter SI.AbsolutePressure p_high = CV.from_bar(250);
      parameter SI.MassFlowRate m_des(fixed = false);
      //Off design calculations
      parameter Integer N_points = 10;
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false, each h.start = 550000);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 550000);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false);
      parameter SI.Power[N_points + 1] W_comp(each fixed = false);
    initial equation
//On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p * PR, state_in[1].h + (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[1].p * PR, s_in[1])) - state_in[1].h) / eta_design);
      W_comp[1] = m_des * (state_out[1].h - state_in[1].h);
//Off-design
      for i in 2:N_points + 1 loop
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_out[i] = MedPB.setState_phX(state_in[i].p * PR, state_in[i].h + (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[i].p * PR, s_in[i])) - state_in[i].h) / eta_design);
        W_comp[i] = m_des * (state_out[i].h - state_in[i].h);
      end for;
    equation

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