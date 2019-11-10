within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign;

 model Exchanger "0D model of a heat exchanger between the HTF and the CO2"
  extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    outer Modelica.Fluid.System system;
    parameter Real pinch = 10;
    parameter Integer N_exch=5;
    parameter Modelica.SIunits.ThermodynamicTemperature TIT_CO2 = 715+273.15;
    parameter Modelica.SIunits.ThermodynamicTemperature T_amb = 40 + 273.15;
    Modelica.Fluid.Interfaces.FluidPort_a HTF_port_a(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {60, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_a CO2_port_a(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {-60, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b HTF_port_b(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b CO2_port_b(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {58, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MedPB.ThermodynamicState[N_exch] state_CO2;
    MedRec.ThermodynamicState[N_exch] state_HTF;
    Modelica.SIunits.SpecificEnthalpy[N_exch] h_CO2(each start = 10^6);
    Modelica.SIunits.SpecificEnthalpy[N_exch] h_HTF;
    SolarTherm.Types.Conductance UA_HX "Conductance of the overall HX";
    SolarTherm.Types.Conductance[N_exch-1] UA_HX_dis "Conductance of the overall HX";
    Real[N_exch] deltaT "Temperature difference in the heat exchangers";
    
    Modelica.SIunits.HeatFlowRate Q_HX_dis;
    Modelica.SIunits.HeatFlowRate Q_HX;
    Modelica.SIunits.ThermodynamicTemperature T_HTF_out;
    SolarTherm.Types.SpecificExergy ex_d;
  protected
  equation
    for i in 1:N_exch loop
      deltaT[i] = MedRec.temperature(state_HTF[i]) - MedPB.temperature(state_CO2[i]);
      state_CO2[i] = MedPB.setState_phX(CO2_port_a.p, h_CO2[i]);
      state_HTF[i] = MedRec.setState_phX(HTF_port_a.p, h_HTF[i]);
    end for;
    h_CO2[N_exch] = CO2_port_b.h_outflow;
    h_HTF[N_exch] = inStream(HTF_port_a.h_outflow);
    h_CO2[1] = inStream(CO2_port_a.h_outflow);
    h_HTF[1] = HTF_port_b.h_outflow;
    TIT_CO2 = MedPB.temperature(state_CO2[N_exch]);
    T_HTF_out = MedRec.temperature(state_HTF[1]);
    
    for i in 1:(N_exch-1) loop
      Q_HX_dis = CO2_port_a.m_flow * (h_CO2[i+1] - h_CO2[i]);
      Q_HX_dis = HTF_port_a.m_flow * (h_HTF[i+1] - h_HTF[i]);
      Q_HX_dis = UA_HX_dis[i] * (deltaT[i] + deltaT[i+1]) / 2;
    end for;
    UA_HX=sum(UA_HX_dis);
    Q_HX=Q_HX_dis*(N_exch-1);
    
    CO2_port_a.p = CO2_port_b.p;
    HTF_port_a.p = HTF_port_b.p;
    CO2_port_a.h_outflow = inStream(CO2_port_b.h_outflow);
    HTF_port_a.h_outflow = inStream(HTF_port_b.h_outflow);
//Closure equation on the mass flow
    CO2_port_a.m_flow = HTF_port_a.m_flow * MedRec.specificHeatCapacityCp(state_HTF[N_exch]) / 1270;
//CO2_port_a.m_flow+CO2_port_b.m_flow=0;
    HTF_port_a.m_flow + HTF_port_b.m_flow = 0;
    ex_d = HTF_port_a.m_flow * (MedRec.specificEnthalpy(state_HTF[N_exch]) - T_amb * MedRec.specificEntropy(state_HTF[N_exch])) + HTF_port_b.m_flow * (MedRec.specificEnthalpy(state_HTF[1]) - T_amb * MedRec.specificEntropy(state_HTF[1])) + CO2_port_b.m_flow * (MedPB.specificEnthalpy(state_CO2[N_exch]) - T_amb * MedPB.specificEntropy(state_CO2[N_exch])) + CO2_port_a.m_flow * (MedPB.specificEnthalpy(state_CO2[1]) - T_amb * MedPB.specificEntropy(state_CO2[1]));
    annotation(
      Diagram(graphics = {Rectangle(origin = {1, 4}, extent = {{-57, 40}, {57, -40}}), Text(origin = {-1, 8}, extent = {{-47, 16}, {47, -16}}, textString = "Exchanger")}),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      Icon(graphics = {Rectangle(origin = {-8, -7}, extent = {{-62, 39}, {62, -39}}), Text(origin = {-6, -13}, extent = {{-48, 21}, {48, -21}}, textString = "Exchanger")}));
    annotation(
      Documentation(info = "<html>
  		<p>The exchanger is divided in sub-HX, which is not necessary to be very accurate as the Cp of CO2 is almost constant in this (p,T) region.</p>
  <p>Closure equations are based on the equality of the energy exchanged and by imposing an outlet temperature of the CO2: the TIT.</p>
  <p>Parameter to integrate in the off-design PB are the UA_HX_dis.</p>
  		</html>"));
  end Exchanger;