within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex;

 model Turbine "OD model of a turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter SI.Efficiency eta_design = 0.9 "isentropic efficiency of the turbine";
      parameter SI.Efficiency PR = 3 "Pressure ratio";
      parameter SI.AngularVelocity N_shaft = 3358;
      parameter SI.Power P_nom = 100 * 10 ^ 6;
      parameter SI.AbsolutePressure p_high = 250 * 10 ^ 5;
      //parameters for on-design
      parameter Modelica.SIunits.Area A_nozzle(fixed = false);
      parameter SI.Diameter diam_turb(fixed = false);
      parameter SI.Velocity tipSpeed(fixed = false);
      //parameters for off and on-design
      parameter Integer N_points = 10;
      parameter SI.MassFlowRate[N_points+1] m_flow(each fixed = false, each start = P_nom / 10 ^ 5);
      parameter SI.Velocity[N_points+1] C_spouting(each fixed = false, each start = 600);
      parameter SI.SpecificEntropy[N_points+1] s_in(each fixed = false, each start = 2900) "inlet entropy to the turbine";
      parameter MedPB.ThermodynamicState[N_points+1] state_in(each p.fixed = false, each h.fixed = false, each p.start = p_high, each h.start = 1.2 * 10 ^ 6);
      parameter MedPB.ThermodynamicState[N_points+1] state_isen(each p.fixed = false, each h.fixed = false, each p.start = p_high / PR);
      parameter SI.AbsolutePressure[N_points] p_in(each fixed = false, each start = p_high);
      parameter SI.ThermodynamicTemperature[N_points] T_in(each fixed = false);
      parameter SI.AbsolutePressure[N_points] p_out(each fixed = false, each start = p_high / PR);
      parameter SI.Power[N_points+1] W_turb(each fixed = false);
      parameter SI.Efficiency[N_points] eta_turb(each fixed = false, each start = eta_design);
      parameter MedPB.ThermodynamicState[N_points+1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
    initial equation
//Part I. On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_isen[1] = MedPB.setState_psX(p_high / PR, s_in[1]);
      state_out[1] = MedPB.setState_phX(p_high / PR, state_in[1].h + eta_design * (state_isen[1].h - state_in[1].h));
      W_turb[1] = m_flow[1] * (state_out[1].h - state_in[1].h);
      C_spouting[1] ^ 2 = 2 * (state_in[1].h - state_isen[1].h);
      m_flow[1] = C_spouting[1] * A_nozzle * MedPB.density(state_out[1]);
      tipSpeed = N_shaft * diam_turb / 2;
      tipSpeed / C_spouting[1] = 0.707;
//Part II. Off-design
      for i in 2:N_points+1 loop
        state_in[i] = MedPB.setState_pTX(p_in[i - 1], T_in[i - 1]);
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_isen[i] = MedPB.setState_psX(p_out[i - 1], s_in[i]);
        state_out[i] = MedPB.setState_phX(p_out[i - 1], state_in[i].h + eta_turb[i - 1] * (state_isen[i].h - state_in[i].h));
        C_spouting[i] ^ 2 = 2 * (state_in[i].h - state_isen[i].h);
        m_flow[i] = C_spouting[i] * A_nozzle * MedPB.density(state_out[i]);
        eta_turb[i - 1] = eta_design * 2 * (tipSpeed / C_spouting[i]) * sqrt(1 - (tipSpeed / C_spouting[i]) ^ 2);

        W_turb[i] = m_flow[i] * (state_out[i].h - state_in[i].h);
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