within SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign;

model simpleRecupPB
    import FI = SolarTherm.Models.Analysis.Finances;
    import SI = Modelica.SIunits;
    extends SolarTherm.Media.CO2.PropCO2;
    extends Icons.PowerBlock;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    input SI.ThermodynamicTemperature T_amb;
      Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation (Placement(
          transformation(extent={{-54,22},{-34,42}}),  iconTransformation(extent={{-48,30},
              {-40,38}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation (Placement(
          transformation(extent={{-74,-60},{-54,-40}}),iconTransformation(extent={{-62,-48},
              {-54,-40}})));
    
    Modelica.Blocks.Interfaces.RealOutput W_net(
      quantity="Power",
      unit="W",
      displayUnit="MW") "Net electric power output" annotation (Placement(
          visible = true,transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
    // PB parameters
    
    parameter Boolean external_parasities = false "= true enable parasities as an input";
    parameter Real nu_min=0.25 "Minimum turbine operation" ;
    Modelica.Blocks.Interfaces.RealInput parasities if external_parasities annotation (Placement(
          transformation(extent={{-12,-12},{12,12}},
          rotation=-90,
          origin={1.77636e-015,80}),                  iconTransformation(
          extent={{-6,-6},{6,6}},
          rotation=-90,
          origin={20,60})));
   
    //Cycle parameters
    parameter SI.AbsolutePressure p_high = 250 * 10 ^ 5 "high pressure of the cycle";
    parameter SI.ThermodynamicTemperature T_high = 800 +273.15 "inlet temperature of the particles";
    parameter SI.ThermodynamicTemperature TIT_CO2 = 715 +273.15 "inlet temperature of the turbine";
    parameter SI.ThermodynamicTemperature T_amb_des = 40 + 273.15 "ambiant temperature";
    parameter Real PR = 3.125 "Pressure ratio";
    parameter SI.Power P_gro = 100*10^6 "first guess of power outlet";
    parameter SI.Power P_nom (fixed=false) "Electrical power at design point";
    parameter SI.MassFlowRate m_HTF_des = 766 "Mass flow rate at design point";
    // Compressor parameters
    parameter SI.Efficiency eta_comp = 0.89 "Maximal isentropic efficiency of the compressors";
    parameter SI.AngularVelocity[4] choiceN = {75000,30000,10000,3600}*0.10471975512 ;
    parameter SI.AngularVelocity N_shaft=(choiceN[integer(Modelica.Math.log(P_gro/10^6)/Modelica.Math.log(10))+2]);
    
    //Turbine parameters
    parameter SI.Efficiency eta_turb = 0.93 "Maximal isentropic efficiency of the turbine";
    
    //Heat recuperator parameters
    parameter Integer N_q = 15;
    
    //Exchanger parameter
    parameter Integer N_exch=5;
    
    //Cooler parameters
    parameter SI.ThermodynamicTemperature T_low = 45 + 273.15 "Outlet temperature of the cooler";
    parameter Integer N_cooler=15;
    
    //Financial analysis
    parameter FI.Money C_HTR (fixed=false) "cost of the heat recuperator";
    parameter FI.Money C_turbine (fixed=false) "cost of the turbine";
    parameter FI.Money C_compressor (fixed=false) "cost of the compressor";
    parameter FI.Money C_exchanger (fixed=false) "cost of the exchanger";
    parameter FI.Money C_generator (fixed=false) "cost of the generator";
    parameter FI.Money C_cooler (fixed=false) "cost of the cooler";
    parameter FI.Money C_PB (fixed=false) "Overall cost of the power block";
    parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
    
    //Results
      SI.Efficiency eta_cycle;

      Boolean m_sup "Disconnect the production of electricity when the outlet pressure of the turbine is close to the critical pressure";
      
    //Components instanciation
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.HeatRecuperatorDTAve HTR(N_q = N_q, P_nom_des = P_nom, pinchRecuperator = 5) annotation(
      Placement(visible = true, transformation(origin = {-38, -16}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.CompressorOnShaft compressor(N_design = N_shaft, eta_design = eta_comp, P_nom_des = P_nom,PR=PR) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Cooler cooler(T_low = T_low, P_nom_des = P_nom,T_amb_des=T_amb_des,N_cool=N_cooler) annotation(
      Placement(visible = true, transformation(origin = {-78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Turbine turbine(PR = PR,  N_shaft = N_shaft, eta_design = eta_turb) annotation(
      Placement(visible = true, transformation(origin = {16, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Exchanger exchanger(P_nom_des = P_nom,T_out_CO2_des=TIT_CO2) annotation(
      Placement(visible = true, transformation(origin = {-4, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = T_high, p_out = 10 ^ 5, redeclare package MedPB = MedRec, use_m_parameter = true, m_flow=turbine.m_des) annotation(
        Placement(visible = true, transformation(origin = {34, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {-36, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      
    
  initial equation
    exchanger.state_HTF_des[N_exch] = MedRec.setState_pTX(10^5,T_high);
    exchanger.m_HTF_des=m_HTF_des;
    P_nom = (-turbine.W_turb_des) - compressor.W_comp_des - cooler.P_cool_des;
// thermodynamic state equalities
    HTR.state_comp_des[1] = compressor.state_out_des;
    exchanger.state_CO2_des[1] = HTR.state_comp_des[N_q];
    turbine.state_in_des = exchanger.state_CO2_des[N_exch];
    turbine.state_out_des = HTR.state_turb_des[N_q];
    cooler.state_des[1] = HTR.state_turb_des[1];
    compressor.state_in_des = cooler.state_des[N_cooler];
  
//mass flow equalities
    turbine.m_des = exchanger.m_CO2_des;
    HTR.m_comp_des = compressor.m_des;
    HTR.m_turb_des = turbine.m_des;
    cooler.m_des = HTR.m_turb_des;
    compressor.m_des = cooler.m_des;
    
// Financial Analysis
    C_HTR= if MedPB.temperature(turbine.state_out_des)>=550+273.15 then 49.45*HTR.UA_HTR^0.7544*(1+0.02141*(MedPB.temperature(turbine.state_out_des)-550-273.15)) else 49.45*HTR.UA_HTR^0.7544;
    C_turbine= if TIT_CO2>= 550+273.15 then 406200*(-turbine.W_turb_des/10^6)^0.8*(1+1.137*10^(-5)*(TIT_CO2-550-273.15)^2) else 406200*(-turbine.W_turb_des/10^6)^0.8;
    C_compressor = 1230000*(compressor.W_comp_des/10^6)^0.3392;
    C_cooler = 32.88*cooler.UA_cooler^0.75;
    C_generator = 108900*(P_nom/10^6)^0.5463;
    C_exchanger = pri_exchanger*exchanger.Q_HX_des/1000;
    C_PB=(C_HTR+C_turbine+C_compressor+C_generator+C_cooler+C_exchanger)*1.05;
    // 1.05 corresponds to inflation from 2017, as correlations are in 2017' dollars.
  equation
    connect(fluid_b, exchanger.HTF_port_b) annotation(
      Line(points = {{-64, -50}, {12, -50}, {12, 10}, {38, 10}, {38, 10}}, color = {0, 127, 255}));
    connect(fluid_a, exchanger.HTF_port_a) annotation(
      Line(points = {{-44, 32}, {54, 32}, {54, 10}, {52, 10}, {52, 10}}, color = {0, 127, 255}));
    connect(src.port_b, fluid_a) annotation(
      Line(points = {{-60, 64}, {-72, 64}, {-72, 30}, {-42, 30}, {-42, 32}, {-44, 32}}, color = {0, 127, 255}));
    connect(sink.port_a, fluid_b) annotation(
      Line(points = {{-42, -78}, {-64, -78}, {-64, -54}, {-64, -54}, {-64, -50}}, color = {0, 127, 255}));
  connect(exchanger.CO2_port_b, turbine.port_a) annotation(
    Line(points = {{2, 6}, {14, 6}, {14, -18}, {12, -18}, {12, -20}}, color = {0, 127, 255}));
  connect(exchanger.CO2_port_a, HTR.from_comp_port_b) annotation(
    Line(points = {{-12, 6}, {-28, 6}, {-28, -8}, {-28, -8}}, color = {0, 127, 255}));
  connect(cooler.port_b, compressor.port_a) annotation(
    Line(points = {{-78, -46}, {-88, -46}, {-88, -6}, {-82, -6}, {-82, -6}}, color = {0, 127, 255}));
  connect(HTR.from_turb_port_b, cooler.port_a) annotation(
    Line(points = {{-49.52, -25.28}, {-49.52, -25.28}, {-49.52, -61.28}, {-77.52, -61.28}, {-77.52, -61.28}}, color = {0, 127, 255}));
  connect(compressor.port_b, HTR.from_comp_port_a) annotation(
    Line(points = {{-71, -17.4}, {-59, -17.4}, {-59, -9.4}, {-51, -9.4}, {-51, -9.4}}, color = {0, 127, 255}));
  connect(turbine.port_b, HTR.from_turb_port_a) annotation(
    Line(points = {{22, -30}, {-28, -30}, {-28, -25}}, color = {0, 127, 255}));
  connect(cooler.T_amb, T_amb);
  connect(m_sup, exchanger.m_sup);
  
  m_sup = exchanger.HTF_port_a.m_flow>= nu_min*exchanger.m_HTF_des;
  //Closure equation.
  if m_sup then 
    turbine.p_out=compressor.p_out/PR;
  else
    exchanger.CO2_port_a.m_flow=exchanger.m_CO2_des;
  end if;
  
  eta_cycle = W_net / exchanger.Q_HX;
  W_net = if m_sup then (-turbine.W_turb) - compressor.W_comp else 0;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_NLS_V,LOG_STATS", outputFormat = "mat", s = "dassl"));
    annotation(
      Documentation(info = "<html>
  		<p>Direct-design model of a simple sCO2 cycle with recuperation. On-design is performed in the initial equation section and off-design in the equation section.</p>
  <p>The mass flow in off-design is imposed by the turbine. The turbomachinery has the same set of variables: N, m_flow, p_out, h_out, with closure equation for 2 (performance maps). </p>
  <p> Therefore, one has to provide values for the two others. 
  <ul>
  <li>Compressors: N, m_flow provided. Outlet: p_out, h_out.</li>
  <li>Turbine: N, p_out provided. Outlet: m_flow, h_out. </li></ul>
   </p>
   <p> p_out,turb=p_out,comp/PR : this closure equation appeared to be the most logical, as one can interfere with a power block at full-load with the outlet pressure of the turbine, and that keeping the PR constant often augments the efficiency. It obviously has to be discussed </p>
  
  <p> The currency is 2019$. The price of the exchanger is taken at 150$/kW_th because it is defined as the objective to reach for next-Gen CSP with particles</p>
  
  <p>N. T. Weiland, B. W. Lance, et S. R. Pidaparti, « SCO2 power cycle components cost correlations from DOE data spanning multiple scales and application », p. 17.</p>
  <p> Available at https://www.netl.doe.gov/projects/files/sCO2PowerCycleComponentCostCorrelationsfromDOEDataSpanningMultipleScalesandApplications_061819.pdf </p>
  <p>&nbsp;</p>
  		</html>"));
  end simpleRecupPB;