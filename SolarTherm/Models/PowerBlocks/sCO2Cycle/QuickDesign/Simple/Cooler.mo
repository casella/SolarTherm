within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple;

 model Cooler
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      input Boolean m_sup;
      input SI.ThermodynamicTemperature T_amb "Ambiant temperature in Kelvin";
      SI.Power P_cooling "Cooling power necessary to cool down the fluid";
      parameter SI.ThermodynamicTemperature T_amb_des = 40 + 273.15 "Ambiant temperature in Kelvin at design point";
      parameter SI.ThermodynamicTemperature T_low = 55 + 273.15;
      parameter SI.Power P_nom_des = 164000;
      parameter Integer N_cooler = 15;
      //On-design only
      parameter MedPB.ThermodynamicState[N_cooler] state_cooler(each p.fixed = false, each h.fixed = false, each h.start = 4.5 * 10 ^ 5);
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SI.Power P_cool_des(fixed = false) "on-design power necessary to run the fans";
      parameter SI.HeatFlowRate Q_dis_cooler(fixed = false);
      parameter SolarTherm.Types.Conductance[N_cooler - 1] UA_cooler_dis(each fixed = false);
      parameter SolarTherm.Types.Conductance UA_cooler(fixed = false) "Conductance of the cooler in W/K";
      parameter MedPB.ThermodynamicState state_out(p.fixed = false, h.fixed = false);
    initial equation
      for k in 1:N_cooler - 1 loop
        Q_dis_cooler = state_cooler[k + 1].h - state_cooler[k].h;
        m_des * Q_dis_cooler = -UA_cooler_dis[k] * (MedPB.temperature(state_cooler[k + 1]) - T_amb_des + MedPB.temperature(state_cooler[k]) - T_amb_des) / 2;
        if k > 1 then
          state_cooler[k].p = state_cooler[1].p;
        end if;
      end for;
      state_cooler[N_cooler] = state_out;
      UA_cooler = sum(UA_cooler_dis);
      P_cool_des * (T_low - T_amb_des) / (-m_des * Q_dis_cooler * (N_cooler - 1)) = 1.49 * 10 ^ 6 * (35.7 - 30) / (136.6 * 10 ^ 6);
    equation
      P_cooling = P_cool_des * ((T_low - T_amb_des) / (max(T_amb + 5, T_low) - T_amb)) ^ (3 / 0.805);
      annotation(
        Documentation(info = "<html>
  		<p>The cooler is thought to be a dry-air cooling device. The outlet temperature of the CO2 is imposed as max(T_low_cycle,T_amb+3). The variation of the ambiant temperature is taken into account in the estimation of the electricity demand for the fans, such as: P_cooling*deltaT/Q_cooler is a constant, deltaT being the average of the temperature of the CO2 and the ambiant, and Q_cooler the energy to withdraw.</p>
  		</html>"));
      annotation(
        Icon(graphics = {Rectangle(origin = {2, 1}, extent = {{-58, 65}, {58, -65}}), Text(origin = {0, -1}, extent = {{-40, -15}, {40, 15}}, textString = "COOLER")}),
  Diagram(graphics = {Rectangle(origin = {-4, 7}, extent = {{-64, 67}, {64, -67}}), Text(origin = {5, 14}, extent = {{-41, -12}, {41, 12}}, textString = "COOLER")}));
    end Cooler;