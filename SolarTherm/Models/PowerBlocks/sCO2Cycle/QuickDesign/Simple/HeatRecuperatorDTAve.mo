within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple;

 model HeatRecuperatorDTAve "The heat recuperator is subdivised in N_q segments in order to accurately represent the CO2 properties variation."
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter Integer N_q = 15 "Number of subdivision of the HX";
      parameter Real ratio_m_des = 1 "ratio of m_comp_des/m_turb_des; we suppose m_turb_des=1, and then scale-up";
      parameter Real pinchRecuperator = 5 "pinch of the recuperator. Imposed as a closing equation for on-design";
      parameter SI.Power P_nom = 100 * 10 ^ 6 "on-design power, useful for start values";
      //On-design
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SolarTherm.Types.Conductance[N_q - 1] UA_dis(each fixed = false, each start = P_nom * 0.4 / N_q);
      parameter SI.HeatFlowRate Q_dis(fixed = false, start = 2 * 10 ^ 4);
      //Off-design
      parameter Integer N_points = 10;
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_turb(each p.fixed = false, each h.fixed = false, each h.start = 7 * 10 ^ 5);
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_comp(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
      parameter SI.TemperatureDifference[N_points + 1, N_q] deltaT(each start = 20, each fixed = false);
    initial equation
//On-design
      for k in 1:N_q - 1 loop
        Q_dis = state_turb[1, k + 1].h - state_turb[1, k].h;
        Q_dis = ratio_m_des * (state_comp[1, k + 1].h - state_comp[1, k].h);
        deltaT[1, k] = MedPB.temperature(state_turb[1, k]) - MedPB.temperature(state_comp[1, k]);
        m_des * Q_dis = UA_dis[k] * (deltaT[1, k] + deltaT[1, k + 1]) / 2;
        state_comp[1, k + 1].p = state_comp[1, 1].p;
        state_turb[1, k].p = state_turb[1, N_q].p;
      end for;
      deltaT[1, N_q] = MedPB.temperature(state_turb[1, N_q]) - MedPB.temperature(state_comp[1, N_q]);
      min(deltaT[1]) = pinchRecuperator;
//Off-design
      for i in 2:N_points + 1 loop
        for k in 1:N_q - 1 loop
          state_turb[i, k + 1].h - state_turb[i, k].h = ratio_m_des * (state_comp[i, k + 1].h - state_comp[i, k].h);
          deltaT[i, k] = MedPB.temperature(state_turb[i, k]) - MedPB.temperature(state_comp[i, k]);
          m_des * (state_turb[i, k + 1].h - state_turb[i, k].h) = UA_dis[k] * (deltaT[i, k] + deltaT[i, k + 1]) / 2;
          state_turb[i, k].p = state_turb[i, N_q].p;
          state_comp[i, k + 1].p = state_comp[i, 1].p;
        end for;
        deltaT[i, N_q] = MedPB.temperature(state_turb[i, N_q]) - MedPB.temperature(state_comp[i, N_q]);
      end for;
      annotation(
        Diagram(graphics = {Rectangle(origin = {1, 7}, extent = {{-61, 31}, {61, -31}}), Text(origin = {5, 1}, extent = {{-53, -17}, {53, 17}}, textString = "RECUPERATOR")}),
  Icon(graphics = {Rectangle(origin = {-3, -9}, extent = {{-65, 33}, {65, -33}}), Text(origin = {-2, -5}, extent = {{-46, -15}, {46, 15}}, textString = "RECUPERATOR")}));
      annotation(
        Documentation(info = "<html>
  		<p>This heat recuperator is a counter-flow HX. Closure equations are based on the equality of m_flow*delta_H for both sides and m_flow*delta_H= UA_i*DTAve_i, DTAve being the average of the temperature difference between the inlet and the outlet of the sub-HX.</p>
  <p>The UA_i must be given as parameters from the on-design analysis.&nbsp;</p>
  		
  		</html>"));
    end HeatRecuperatorDTAve;