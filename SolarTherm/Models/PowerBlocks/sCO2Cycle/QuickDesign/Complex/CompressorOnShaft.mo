within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex;

 model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      parameter SI.Efficiency eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
      parameter SI.AngularVelocity N_design = 3358 "Design rotationnal speed in rad/s";
      parameter Real PR = 3 "pressure ratio chosen";
      //On-design calculations
      parameter SI.AbsolutePressure p_high = CV.from_bar(250);
      parameter SI.MassFlowRate[N_points+1] m_flow(each fixed = false);
      parameter SI.Diameter diam_rotor(fixed = false) "on-design diameter of the rotor";
      parameter SI.Velocity tipSpeed(fixed = false) " tip speed of the rotor";
      parameter Real psi_des(fixed = false) "on-design adimensionned head";
      parameter Real phi_opt = 0.0297035 "optimal adimensionned mass flow";
      //Off design calculations
      parameter Integer N_points = 10;
      parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
      parameter MedPB.ThermodynamicState[N_points+1] state_in(each p.fixed = false, each h.fixed = false, each h.start = 450000, each p.start = p_high / PR);
      parameter MedPB.ThermodynamicState[N_points+1] state_isen(each p.fixed = false, each h.fixed = false, each p.start = p_high, each h.start = 500000);
      parameter MedPB.ThermodynamicState[N_points+1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 550000, each p.start = p_high);
      parameter SI.SpecificEntropy[N_points+1] s_in(each fixed = false, each start = 1500);
      parameter Real[N_points] psi(each fixed = false, each start = psi_des) "off-design adimensionned head";
      parameter Real[N_points] phi(each fixed = false, each start = phi_opt) "off-design adimensionned mass flow";
      parameter SI.AbsolutePressure[N_points] p_out(each fixed = false, each start = p_high);
      parameter SI.Efficiency[N_points] eta_comp(each fixed = false, each start = eta_design);
      parameter SI.Power[N_points+1] W_comp(each fixed = false);
    initial equation
//On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_isen[1] = MedPB.setState_psX(state_in[1].p * PR, s_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p * PR, state_in[1].h + (state_isen[1].h - state_in[1].h) / eta_design);
      W_comp[1] = m_flow[1] * (state_out[1].h - state_in[1].h);
      2 * m_flow[1] = phi_opt * MedPB.density(state_in[1]) * N_design * diam_rotor ^ 3;
      tipSpeed = N_design * diam_rotor / 2;
      psi_des = (state_isen[1].h - state_in[1].h) / tipSpeed ^ 2;
//Off-design
      for i in 2:N_points+1 loop
        
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_isen[i] = MedPB.setState_psX(p_out[i - 1], s_in[i]);
        state_out[i] = MedPB.setState_phX(p_out[i - 1], state_in[i].h + (state_isen[i].h - state_in[i].h) / eta_comp[i - 1]);
        W_comp[i] = m_flow[i] * (state_out[i].h - state_in[i].h);
        phi[i - 1] = m_flow[i] / (MedPB.density(state_in[i]) * tipSpeed * diam_rotor ^ 2);
        psi[i - 1] = (state_isen[i].h - state_in[i].h) / tipSpeed ^ 2;
        psi[i - 1] = (0.04049 + 54.7 * phi[i - 1] - 2505 * phi[i - 1] ^ 2 + 53224 * phi[i - 1] ^ 3 - 498626 * phi[i - 1] ^ 4) * psi_des / 0.46181921979961293;
        eta_comp[i - 1] = eta_design / 0.677837 * ((-0.7069) + 168.6 * phi[i - 1] - 8089 * phi[i - 1] ^ 2 + 182725 * phi[i - 1] ^ 3 - 1.638 * 10 ^ 6 * phi[i - 1] ^ 4);

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