within SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign;

  model HeatRecuperatorDTAve "The heat recuperator is subdivised in N_q segments in order to accurately represent the CO2 properties variation."
    extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    import SI = Modelica.SIunits;
    parameter Integer N_q = 15 "Number of subdivision of the HX";
    parameter Real ratio_m_des =1 "ratio of m_comp_des/m_turb_des; we suppose m_turb_des=1, and then scale-up";
    parameter Real pinchRecuperator = 5 "pinch of the recuperator. Imposed as a closing equation for on-design";
      Modelica.Fluid.Interfaces.FluidPort_a from_comp_port_a(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_a from_turb_port_a(redeclare package Medium = MedPB, m_flow.start=P_nom_des/10^5) annotation(
        Placement(visible = true, transformation(origin = {60, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b from_comp_port_b(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b from_turb_port_b(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {-60, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
      Real[N_q] deltaT(start = {150 for i in 1:N_q});
    SI.HeatFlowRate Q_HX;
    parameter SolarTherm.Types.Conductance UA_HTR(fixed = false) "Conductance of the overall HX";
    parameter SI.Power P_nom_des = 10 ^ 5;
    parameter SolarTherm.Types.Conductance UA_dis[N_q - 1](each fixed = false, each start = 0.04 * P_nom_des) "on-design conductance of the heat recuperator";
    
    protected
    parameter MedPB.ThermodynamicState[N_q] state_turb_des(each p.fixed = false, each h.fixed = false,  h.start = {530000 + i / N_q * 500000 for i in 1:N_q}), state_comp_des(each p.fixed = false, each h.fixed = false, h.start = {500000 + i / N_q * 500000 for i in 1:N_q});
  
    parameter SI.TemperatureDifference[N_q] deltaT_des(each fixed = false, each start = 25);
    parameter SI.MassFlowRate m_comp_des(fixed = false, start = P_nom_des / 10 ^ 5), m_turb_des(fixed = false, start = P_nom_des / 10 ^ 5) "on-design mass flow from the compressor, turbine";
    parameter SI.HeatFlowRate Q_HX_des(fixed = false, start=75000);
    parameter SI.HeatFlowRate Q_dis_des(fixed = false, start=18600);
  
  protected
   MedPB.ThermodynamicState[N_q] state_from_turb (h.start = {530000 + i / N_q * 500000 for i in 1:N_q}), state_from_comp (h.start = {500000 + i / N_q * 500000 for i in 1:N_q});
  initial equation
    for i in 2:N_q loop
      deltaT_des[i] = MedPB.temperature(state_turb_des[i]) - MedPB.temperature(state_comp_des[i]);
      state_comp_des[i].p = state_comp_des[1].p;
    end for;
    deltaT_des[1] = MedPB.temperature(state_turb_des[1]) - MedPB.temperature(state_comp_des[1]);
    min(deltaT_des) = pinchRecuperator;
    
    Q_HX_des = m_turb_des*Q_dis_des * (N_q - 1);
    UA_HTR = sum(UA_dis);
    for i in 1:N_q - 1 loop
      state_turb_des[i].p = state_turb_des[N_q].p;
      Q_dis_des = ratio_m_des * (state_comp_des[i + 1].h - state_comp_des[i].h);
      m_turb_des*Q_dis_des = UA_dis[i] * (deltaT_des[i + 1] + deltaT_des[i]) / 2;
      Q_dis_des = (state_turb_des[i + 1].h - state_turb_des[i].h);
    end for;
  equation
  for i in 1:N_q loop
    deltaT[i] = MedPB.temperature(state_from_turb[i])-MedPB.temperature(state_from_comp[i]);
    state_from_turb[i].p = from_turb_port_a.p;
    state_from_comp[i].p = from_comp_port_a.p;
  end for;
  
  state_from_comp[1].h=inStream(from_comp_port_a.h_outflow);
  state_from_turb[N_q].h=inStream(from_turb_port_a.h_outflow);
  from_turb_port_b.h_outflow = state_from_turb[1].h;
  from_comp_port_b.h_outflow = state_from_comp[N_q].h;
  Q_HX=from_turb_port_a.m_flow*(state_from_turb[N_q].h-state_from_turb[1].h);
  for i in 2:N_q loop
    from_turb_port_a.m_flow * (state_from_turb[i].h - state_from_turb[i - 1].h) = from_comp_port_a.m_flow * (state_from_comp[i].h - state_from_comp[i - 1].h);
    from_turb_port_a.m_flow * (state_from_turb[i].h - state_from_turb[i - 1].h) = UA_dis[i - 1] * (abs(from_comp_port_a.m_flow / m_comp_des + from_turb_port_a.m_flow / m_turb_des) ^ 0.8)/ (2 ^ 0.8 )* (deltaT[i-1] + deltaT[i]) / 2;
  end for;
  from_turb_port_b.m_flow + from_turb_port_a.m_flow = 0;
  from_comp_port_b.m_flow + from_comp_port_a.m_flow = 0;
  from_turb_port_b.p = from_turb_port_a.p;
  from_comp_port_b.p = from_comp_port_a.p;
  from_turb_port_a.h_outflow = inStream(from_turb_port_b.h_outflow);
  from_comp_port_a.h_outflow = inStream(from_comp_port_b.h_outflow);
    annotation(
      Documentation(info = "<html>
  		<p>This heat recuperator is a counter-flow HX. Closure equations are based on the equality of m_flow*delta_H for both sides and m_flow*delta_H= UA_i*DTAve_i, DTAve being the average of the temperature difference between the inlet and the outlet of the sub-HX.</p>
  		
  		</html>"));
    annotation(
      Diagram(graphics = {Rectangle(origin = {1, 7}, extent = {{-61, 31}, {61, -31}}), Text(origin = {5, 1}, extent = {{-53, -17}, {53, 17}}, textString = "RECUPERATOR")}),
  Icon(graphics = {Rectangle(origin = {-3, -9}, extent = {{-65, 33}, {65, -33}}), Text(origin = {-2, -5}, extent = {{-46, -15}, {46, 15}}, textString = "RECUPERATOR")}));
  end HeatRecuperatorDTAve;