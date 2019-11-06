within SolarTherm.Models.PowerBlocks.sCO2Cycle;

package DirectDesign
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
      SI.Power P_elec;
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
  //    SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = T_high, p_out = 10 ^ 5, redeclare package MedPB = MedRec, use_m_parameter = true, m_flow=turbine.m_des) annotation(
  //      Placement(visible = true, transformation(origin = {34, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  //    SinkFlow sink annotation(
  //      Placement(visible = true, transformation(origin = {-36, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      
    
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
  //  connect(src.port_b, fluid_a) annotation(
  //    Line(points = {{-60, 64}, {-72, 64}, {-72, 30}, {-42, 30}, {-42, 32}, {-44, 32}}, color = {0, 127, 255}));
  //  connect(sink.port_a, fluid_b) annotation(
  //    Line(points = {{-42, -78}, {-64, -78}, {-64, -54}, {-64, -54}, {-64, -50}}, color = {0, 127, 255}));
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
  
  eta_cycle = P_elec / exchanger.Q_HX;
  P_elec = if m_sup then (-turbine.W_turb) - compressor.W_comp else 0;
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

  model Exchanger
    extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    import SI = Modelica.SIunits;
    parameter SI.ThermodynamicTemperature T_out_CO2_des=715+273.15;
    parameter SI.Power P_nom_des = 164000;
      input Boolean m_sup "when m_sup=false, m_HTF=m_HTF_design and P_elec=0 -> allows switching off the PB";
    parameter Real ratio_m_des = 1 "ratio of m_CO2_des/m_HTF_des at design point";  
    parameter Integer N_exch = 8;
      Modelica.Fluid.Interfaces.FluidPort_a HTF_port_a(redeclare package Medium = MedRec) annotation(
        Placement(visible = true, transformation(origin = {60, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_a CO2_port_a(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {-60, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-72, -56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b HTF_port_b(redeclare package Medium = MedRec) annotation(
        Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b CO2_port_b(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {58, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      MedPB.ThermodynamicState[N_exch] state_CO2, state_HTF;
      SI.SpecificEnthalpy[N_exch] h_CO2(start = {990000 + (i/N_exch)*200000 for i in 1:N_exch}),h_HTF (start={600000 + (i/N_exch)*200000 for i in 1:N_exch});
      Real[N_exch] deltaT "Temperature difference in the heat exchangers";
      SI.HeatFlowRate Q_HX;
      SI.ThermodynamicTemperature T_CO2_out, T_HTF_out;
      //Real deltaT_lm;
      Real deltaTAve;
      SI.MassFlowRate m_HTF_bis (start=P_nom_des/10^5);
      
    protected
    parameter SI.HeatFlowRate Q_HX_des(fixed = false, start= 10^5*2/N_exch);
    parameter SI.MassFlowRate m_CO2_des(fixed = false), m_HTF_des(fixed = false);
    parameter SolarTherm.Types.Conductance UA_HX(fixed = false) "on-design conductance of the overall exchanger";
    parameter SolarTherm.Types.Conductance[N_exch-1] UA_HX_dis(each fixed = false) "on-design conductance of the exchanger";
    
    parameter Real[N_exch] deltaT_des(each fixed = false, each start = 75);
  
    parameter MedPB.ThermodynamicState[N_exch] state_CO2_des(each p.fixed = false, each h.fixed = false, each p.start = 250 * 10 ^ 5, each h.start = 10 ^ 6), state_HTF_des(each p.fixed = false, each h.fixed = false, each p.start = 10 ^ 5, each h.start = 855004);
  initial equation
    for i in 2:(N_exch-1) loop
      
      state_CO2_des[i].p = state_CO2_des[1].p;
    end for;
    state_CO2_des[N_exch]=MedPB.setState_pTX(state_CO2_des[1].p,T_out_CO2_des);
    deltaT_des[N_exch] = MedRec.temperature(state_HTF_des[N_exch]) - MedPB.temperature(state_CO2_des[N_exch]);
    for i in 1:(N_exch-1) loop
      deltaT_des[i] = MedRec.temperature(state_HTF_des[i]) - MedPB.temperature(state_CO2_des[i]);
      state_HTF_des[i].p = state_HTF_des[N_exch].p;
      Q_HX_des = ratio_m_des * (state_CO2_des[i+1].h - state_CO2_des[i].h);
      Q_HX_des =  (state_HTF_des[i+1].h - state_HTF_des[i].h);
      m_HTF_des*Q_HX_des = UA_HX_dis[i] * (deltaT_des[i] + deltaT_des[i+1]) / 2;
    end for;
    UA_HX=sum(UA_HX_dis);
    
    m_CO2_des = ratio_m_des * m_HTF_des;
  equation
  for i in 1:N_exch loop
    deltaT[i] = if m_sup then MedRec.temperature(state_HTF[i]) - MedPB.temperature(state_CO2[i]) else deltaT_des[i];
    state_CO2[i] = MedPB.setState_phX(CO2_port_a.p, h_CO2[i]);
    state_HTF[i] = MedRec.setState_phX(HTF_port_a.p, h_HTF[i]);
  end for;
  T_CO2_out = MedPB.temperature(state_CO2[N_exch]);
  T_HTF_out = MedRec.temperature(state_HTF[1]);
  
  deltaTAve = (deltaT[1] + deltaT[N_exch]) / 2;
  h_CO2[N_exch] = CO2_port_b.h_outflow;
  h_HTF[N_exch] = if m_sup then inStream(HTF_port_a.h_outflow) else state_HTF_des[N_exch].h;
  h_CO2[1] = inStream(CO2_port_a.h_outflow);
  HTF_port_b.h_outflow = if m_sup then h_HTF[1] else inStream(HTF_port_a.h_outflow);
  m_HTF_bis = if m_sup then HTF_port_a.m_flow else m_HTF_des;
  Q_HX = CO2_port_a.m_flow * (h_CO2[N_exch] - h_CO2[1]);
  for i in 1:(N_exch-1) loop 
    m_HTF_bis*(h_HTF[i+1]-h_HTF[i])=CO2_port_a.m_flow*(h_CO2[i+1]-h_CO2[i]);
    CO2_port_a.m_flow*(h_CO2[i+1]-h_CO2[i])=UA_HX_dis[i]* (1 / 2 * abs(m_HTF_bis / m_HTF_des + CO2_port_a.m_flow / m_CO2_des)) ^ 0.8* (deltaT[i] + deltaT[i+1]) / 2;
  end for;
  
  HTF_port_a.h_outflow = inStream(HTF_port_b.h_outflow);
  CO2_port_a.h_outflow = inStream(CO2_port_b.h_outflow);
  //It is necessary to have one equation in a cycle that doesn't imply a circular equality on the mass flow rates
  //CO2_port_b.m_flow + CO2_port_a.m_flow = 0;
  HTF_port_a.m_flow + HTF_port_b.m_flow = 0;
  //CO2_port_a.m_flow = if m_sup then HTF_port_a.m_flow else m_CO2_des * 0.8;
  // Pressure equality
  CO2_port_b.p = CO2_port_a.p;
  HTF_port_a.p = HTF_port_b.p;
    annotation(
      Diagram(graphics = {Rectangle(origin = {1, 4}, extent = {{-57, 40}, {57, -40}}), Text(origin = {-1, 8}, extent = {{-47, 16}, {47, -16}}, textString = "Exchanger")}),
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    annotation(
      Documentation(info = "<html>
  		<p>The exchanger is a heat exchanger between the HTF and the CO2. It is a counterflow HX, based on a TLMD. The conductance UA has to be specified from the on-design.</p>
  <p>The conductance in off-design varies as UA_Off=UA_on*(m_flow/m_design)^0.8.&nbsp;<span >The average between the two mass flows is taken.</span></p>
  <p>A.T. Louis et T. Neises, analysis and optimization for Off-design performance of the recompression s-CO2 cycles for high temperature CSP applications, in The 5th International Symposium-Supercritical CO2 Power Cycles, 2016</p>
  <p>&nbsp;</p>
  		</html>"));
  end Exchanger;

  model Cooler
    extends SolarTherm.Media.CO2.PropCO2;
    import SI = Modelica.SIunits;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
  
    input SI.ThermodynamicTemperature T_amb "Ambiant temperature in Kelvin";
    parameter SI.ThermodynamicTemperature T_amb_des = 40 + 273.15 "Ambiant temperature in Kelvin at design point";
    parameter SolarTherm.Types.Conductance UA_cooler(fixed = false) "Conductance of the cooler in W/K";
    parameter SI.ThermodynamicTemperature T_low = 55 + 273.15;
    parameter SI.Power P_nom_des = 164000;
    parameter Integer N_cool = 15;
    parameter Real deltaT_design = 15 "Difference between ambient and outlet CO2 temperature";
    
    MedPB.ThermodynamicState state_a (h.start=500000)"Thermodynamic State at the entrance";
    MedPB.ThermodynamicState state_b "Thermodynamic State at the outlet";
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB, m_flow.start=P_nom_des/10^5) annotation(
        Placement(visible = true, transformation(origin = {0, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {2, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
        Placement(visible = true, transformation(origin = {-2, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real[2] deltaT;
    SI.HeatFlowRate Q_cooler;
    SI.Power P_cooling;
    
  
    
    parameter SolarTherm.Types.Conductance UA_dis[N_cool - 1](each fixed = false) "Conductance of the cooler per sub-HX";
    parameter SI.Power P_cool_des(fixed = false, start = 0.01 * P_nom_des) "power necessary to cool down at design point";
    protected
    parameter SI.MassFlowRate m_des(fixed = false);
    parameter SI.HeatFlowRate Q_cooler_des(fixed = false, start = 10 ^ 6);
    parameter MedPB.ThermodynamicState[N_cool] state_des(each p.fixed = false, each h.fixed = false) "Thermodynamic State at the i-th position";
  
    //parameter Real[N_cool-1] deltaT_lm_des (each fixed=false)"logarithmic temperature difference";
    parameter Real[N_cool] deltaT_des(each fixed = false, each start = 25) "difference with the ambiant air at the inlet and outlet";
    parameter SI.HeatFlowRate Q_dis_des(fixed = false, start = 10 ^ 5) "Heat flow rate dispatched per sub-HX in the cooler";
    parameter SI.ThermodynamicTemperature[N_cool] T_CO2_des(each fixed = false, each start = 273.15 + 75);
  initial equation
    for i in 2:N_cool loop
      state_des[i] = MedPB.setState_pTX(state_des[1].p, T_CO2_des[i]);
      deltaT_des[i] = T_CO2_des[i] - T_amb_des;
    end for;
    T_CO2_des[N_cool] = T_low;
    deltaT_des[1] = T_CO2_des[1] - T_amb_des;
    T_CO2_des[1]=MedPB.temperature(state_des[1]);
    for i in 1:N_cool - 1 loop
      Q_dis_des = (state_des[i+1].h - state_des[i].h);
      m_des * Q_dis_des = -UA_dis[i] * (deltaT_des[i] + deltaT_des[i + 1]) / 2;
    end for;
    UA_cooler = sum(UA_dis);
    Q_cooler_des = (N_cool - 1) * Q_dis_des*m_des;
    P_cool_des * deltaT_des[N_cool]/(-Q_cooler_des)= 1.49*10^6*(35.7-30)/(136.6*10^6);
  
  equation
  deltaT = {MedPB.temperature(state_a) - T_amb, MedPB.temperature(state_b) - T_amb};
  state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
  state_b = MedPB.setState_pTX(port_a.p, max(T_amb + 3, T_low));
  P_cooling = P_cool_des* (deltaT_design / deltaT[2]) ^ (3 / 0.805)*(Q_cooler/Q_cooler_des);
  
  Q_cooler = port_a.m_flow * (state_b.h - state_a.h);
  port_a.m_flow + port_b.m_flow = 0;
  port_a.p = port_b.p;
  port_b.h_outflow = state_b.h;
  port_a.h_outflow = inStream(port_b.h_outflow);
    annotation(
      Icon(graphics = {Rectangle(origin = {2, 1}, extent = {{-58, 65}, {58, -65}}), Text(origin = {0, -1}, extent = {{-40, -15}, {40, 15}}, textString = "COOLER")}),
  Diagram(graphics = {Rectangle(origin = {-4, 7}, extent = {{-64, 67}, {64, -67}}), Text(origin = {5, 14}, extent = {{-41, -12}, {41, 12}}, textString = "COOLER")}));
    annotation(
      Documentation(info = "<html>
  		<p>The cooler is thought to be a dry-air cooling device. The outlet temperature of the CO2 is imposed as max(T_low_cycle,T_amb+5). The variation of the ambiant temperature is taken into account in the estimation of the electricity demand for the fans, such as: P_cooling*deltaT/Q_cooler is a constant, deltaT being the average of the temperature of the CO2 and the ambiant, and Q_cooler the energy to withdraw.</p>
  		</html>"));
  end Cooler;

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

  model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
  extends SolarTherm.Media.CO2.PropCO2;
    import SI = Modelica.SIunits;
    import CV = Modelica.SIunits.Conversions;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    parameter Real eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
    parameter SI.Diameter diam_rotor(fixed = false) "on-design diameter of the rotor";
    parameter SI.AngularVelocity N_design = 40000 * 0.104 "Design rotationnal speed in rad/s";
    parameter Real psi_des(fixed = false) "on-design adimensionned head";
    parameter Real PR = 3 "pressure ratio chosen";
    parameter SI.Power P_nom_des "nominal power at design point";
    parameter SI.ThermodynamicTemperature T_in_des = CV.from_degC(45) "chosen inlet temperature of the compressor at design point";
    parameter SI.AbsolutePressure p_high_des = CV.from_bar(250);
    parameter Real phi_opt = 0.0297035 "optimal adimensionned mass flow";
    parameter SI.AngularVelocity tipSpeed(fixed = false, start = 500) "tip Speed of the rotor";
    SI.AbsolutePressure p_out(start = p_high_des) "outlet pressure";
    MedPB.ThermodynamicState state_a "thermodynamic state at the entrance";
    MedPB.ThermodynamicState state_isen "thermodynamic state at the end of the isentropic compression";
    MedPB.ThermodynamicState state_b "thermodynamic state at the end of the real compresssion";
    SI.Power W_comp "power used for compression";
    SI.SpecificEntropy s_entrance "entropy at the entrance of the compressor";
    Real phi "adimensionned mass flow rate";
    Real psi "adimensionned head";
    SI.Density d_entrance(start = 267) "density at the inlet";
    SI.Efficiency eta_comp(start = eta_design) "isentropic efficiency of the compressor";
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {42, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
    protected
    parameter MedPB.ThermodynamicState state_in_des(p.fixed = false, h.fixed = false, h.start = 550000),state_isen_des(p.fixed = false, h.fixed = false),state_out_des(p.fixed = false, h.fixed = false);
  
    
    parameter SI.Power W_comp_des(fixed = false);
    parameter SI.MassFlowRate m_des(fixed = false, start = 3) "design mass flow rate in kg/s";
  initial equation
    2 * m_des = phi_opt * MedPB.density(state_in_des) * N_design * diam_rotor ^ 3;
    psi_des = (state_isen_des.h - state_in_des.h) / tipSpeed ^ 2;
    W_comp_des = m_des * (state_out_des.h - state_in_des.h);
  
    state_isen_des = MedPB.setState_psX(p_high_des, MedPB.specificEntropy(state_in_des));
    tipSpeed = diam_rotor * N_design / 2;
    state_out_des = MedPB.setState_phX(p_high_des,state_in_des.h + (state_isen_des.h - state_in_des.h) / eta_design);
  equation
    state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
    s_entrance = MedPB.specificEntropy(state_a);
    state_isen = MedPB.setState_psX(p_out, s_entrance);
    state_b = MedPB.setState_phX(p_out, state_a.h + (state_isen.h - state_a.h) / eta_comp);
    port_b.p = state_b.p;
    port_b.h_outflow = state_b.h;
    W_comp = port_a.m_flow * (state_b.h - state_a.h);
    port_a.h_outflow = 0;
    d_entrance = MedPB.density(state_a);
    port_a.m_flow + port_b.m_flow = 0;
    phi = port_a.m_flow / (d_entrance * tipSpeed * diam_rotor ^ 2);
    psi = (state_isen.h - state_a.h) / tipSpeed ^ 2;
    psi = (0.04049 + 54.7 * phi - 2505 * phi ^ 2 + 53224 * phi ^ 3 - 498626 * phi ^ 4) * psi_des / 0.46181921979961293;
    eta_comp = eta_design / 0.677837 * ((-0.7069) + 168.6 * phi - 8089 * phi ^ 2 + 182725 * phi ^ 3 - 1.638 * 10 ^ 6 * phi ^ 4);
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

  model CompressorToDiscuss
    extends SolarTherm.Media.CO2.PropCO2;
    import SI= Modelica.SIunits;
    import CV = Modelica.SIunits.Conversions;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {-60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
      Placement(visible = true, transformation(origin = {42, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter Real eta_design = 0.87 "Maximal isentropic efficiency of the compressor";
    parameter Real PR = 2.313 "Pressure ratio";
    parameter SI.AbsolutePressure p_high_des = 250 * 10 ^ 5 "outlet pressure";
    parameter Real phi_opt = 0.0297035 "optimal adimensionned mass flow";
    
    MedPB.ThermodynamicState state_a "thermodynamic state at the entrance";
    MedPB.ThermodynamicState state_isen "thermodynamic state at the end of the isentropic compression";
    MedPB.ThermodynamicState state_b "thermodynamic state at the end of the real compresssion";
    SI.Power W_comp "power used for compression";
    SI.AbsolutePressure p_out (start=p_out_des);
    SI.SpecificEntropy s_entrance "exergy at the entrance of the compressor";
    
    Real phi;
    Real psi;
    SI.Velocity tipSpeed (start=tipSpeed_des);
    
    parameter SI.AngularVelocity N_design (fixed=false);
    parameter SI.Velocity tipSpeed_des (fixed=false);
    SI.AngularVelocity N_compressor(start = N_design);
    Real d_entrance(start = 344);
    Real eta_comp(start = eta_design);
    parameter SI.SpecificEnthalpy h_out_des(fixed = false), h_in_des(fixed = false) "enthalpy at the outlet of the compressor";
    parameter MedPB.ThermodynamicState state_in_des(p.fixed = false, h.fixed = false);
    parameter MedPB.ThermodynamicState state_isen_des(p.fixed = false, h.fixed = false);
    parameter SI.AbsolutePressure p_in_des(fixed = false), p_out_des(fixed = false);
    parameter SI.Power W_comp_des(fixed = false);
    parameter SI.MassFlowRate m_des(fixed = false, start = 3) "design mass flow rate in kg/s";  
    parameter Real psi_des(fixed = false) "on-design adimensionned head";
    parameter SI.Diameter diam_rotor(fixed = false) "on-design diameter of the rotor";
    Real Nratio;
    Real deltaHi;
    parameter Real deltaHi_des (fixed=false);
  protected
    
  initial equation
    2 * m_des = phi_opt * MedPB.density(state_in_des) * N_design * diam_rotor ^ 3;
    psi_des = (state_isen_des.h - state_in_des.h) / tipSpeed_des ^ 2;
    W_comp_des = m_des * (h_out_des - state_in_des.h);
    state_in_des = MedPB.setState_phX(p_in_des, h_in_des);
    p_out_des = p_high_des;
    tipSpeed_des=diam_rotor*N_design/2;
    state_isen_des = MedPB.setState_psX(p_out_des, MedPB.specificEntropy(state_in_des));
    h_out_des = state_in_des.h + (state_isen_des.h - state_in_des.h) / eta_design;
    deltaHi_des=state_isen_des.h - state_in_des.h;
  equation
    state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
    s_entrance = MedPB.specificEntropy(state_a);
    state_isen = MedPB.setState_psX(p_out, s_entrance);
    state_b = MedPB.setState_phX(p_out, state_a.h + (state_isen.h - state_a.h) / eta_comp);
    port_b.p = state_b.p;
    port_b.h_outflow = state_b.h;
    port_a.h_outflow = 0;
    
    W_comp = port_a.m_flow * (state_b.h - state_a.h);
  port_a.m_flow + port_b.m_flow = 0;


  d_entrance = MedPB.density(state_a);
// Compressor on-design parameters to be calculated
  
  tipSpeed = diam_rotor * N_compressor / 2;
    
  
  phi = port_a.m_flow / (d_entrance * tipSpeed * diam_rotor ^ 2) * (N_compressor / N_design) ^ (1 / 5);
  psi = (state_isen.h - state_a.h)* (( N_compressor/N_design) ^ ((20 * phi) ^ 3)) / tipSpeed ^ 2 ;
Nratio=( N_compressor/N_design);
deltaHi=psi*tipSpeed^2/(( N_compressor/N_design) ^ ((20 * phi) ^ 3));
  psi = (0.04049 + 54.7 * phi - 2505 * phi ^ 2 + 53224 * phi ^ 3 - 498626 * phi ^ 4) * psi_des / 0.46181921979961293;
  eta_comp = eta_design / 0.677837 * (abs((-0.7069) + 168.6 * phi - 8089 * phi ^ 2 + 182725 * phi ^ 3 - 1.638 * 10 ^ 6 * phi ^ 4)) * (N_compressor/N_design) ^ ((20 * phi) ^ 5);

    annotation(
      Diagram(graphics = {Text(origin = {-20, 18}, extent = {{-28, 16}, {42, -46}}, textString = "COMPRESSOR"), Polygon(origin = {-12, 10}, points = {{-42, 40}, {-42, -44}, {42, -70}, {42, 70}, {-42, 40}, {-42, 40}})}, coordinateSystem(initialScale = 0.1)),
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Polygon(origin = {-26, -2}, points = {{-40, 42}, {-42, -48}, {42, -78}, {42, 78}, {-40, 42}}), Text(origin = {-16, 11}, extent = {{-48, -31}, {24, 15}}, textString = "COMPRESSOR")}));
  end CompressorToDiscuss;

  model recompPB
    extends SolarTherm.Media.CO2.PropCO2;
    import SI = Modelica.SIunits;
    import FI = SolarTherm.Models.Analysis.Finances;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    //input SolarTherm.Interfaces.Connectors.WeatherBus wbus;
    extends Icons.PowerBlock;
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
    input SI.ThermodynamicTemperature T_amb;
    //Cycle parameters
    parameter SI.AbsolutePressure p_high = 200 * 10 ^ 5 "high pressure of the cycle";
    parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the particles";
    parameter SI.ThermodynamicTemperature TIT_CO2 = 715 + 273.15 "inlet temperature of the turbine";
    parameter SI.ThermodynamicTemperature T_amb_des = 30 + 273.15 "ambiant temperature";
    parameter Real PR = 2.5 "Pressure ratio";
    parameter SI.Power P_gro = 100*10^6 "first guess of power outlet";
    parameter SI.Power P_nom (fixed=false) "Electrical power at design point";
    parameter SI.MassFlowRate m_HTF_des = 1000 "Mass flow rate at design point";
    parameter Real gamma = 0.28 "Part of the mass flow going to the recompression directly";
    parameter SI.AngularVelocity[4] choiceN = {75000,30000,10000,3600}*0.10471975512 ;
    parameter SI.AngularVelocity N_shaft=(choiceN[integer(Modelica.Math.log(P_gro/10^6)/Modelica.Math.log(10))+2]);
    
    // main Compressor parameters
    parameter SI.Efficiency eta_comp_main = 0.89 "Maximal isentropic efficiency of the compressors";
  
    // reCompressor parameters
    parameter SI.Efficiency eta_comp_re = 0.89 "Maximal isentropic efficiency of the compressors";
  
    //Turbine parameters
    parameter SI.Efficiency eta_turb = 0.93 "Maximal isentropic efficiency of the turbine";
    
    //HTR Heat recuperator parameters
    parameter Integer N_HTR = 15;
  
    //LTR Heat recuperator parameters
    parameter Integer N_LTR = 15;
    parameter Real ratio_m_des=1-gamma;
    //Cooler parameters
    parameter SI.ThermodynamicTemperature T_low = 45 + 273.15 "Outlet temperature of the cooler";
    parameter Integer N_cooler = 15;
    
    //Exchanger parameters
    parameter SI.ThermodynamicTemperature T_HTF_in_des = 800+273.15;
    parameter Integer N_exch = 5;
    //Financial analysis
    parameter FI.Money C_HTR (fixed=false) "cost of the high temperature heat recuperator";
    parameter FI.Money C_LTR (fixed=false) "cost of the low temperature heat recuperator";
    parameter FI.Money C_turbine (fixed=false) "cost of the turbine";
    parameter FI.Money C_mainCompressor (fixed=false) "cost of the main compressor";
    parameter FI.Money C_reCompressor (fixed=false) "cost of the re compressor";
    parameter FI.Money C_exchanger (fixed=false) "cost of the exchanger";
    parameter FI.Money C_generator (fixed=false) "cost of the generator";
    parameter FI.Money C_cooler (fixed=false) "cost of the cooler";
    parameter FI.Money C_PB (fixed=false) "Overall cost of the power block";
    parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
    
    
    //Results
    SI.Efficiency eta_cycle;
    
    SI.Energy E_net(final start=0, fixed=true, displayUnit="MW.h");
    Boolean m_sup "Disconnect the production of electricity when the outlet pressure of the turbine is close to the critical pressure";
    
    //Components instanciation
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.HeatRecuperatorDTAve HTR(N_q = N_HTR,P_nom_des=P_gro, ratio_m_des=1) annotation(
      Placement(visible = true, transformation(origin = {12, -22}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.CompressorOnShaft mainCompressor(eta_design = eta_comp_main, N_design=N_shaft,P_nom_des=P_gro,p_high_des=p_high) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Cooler cooler(T_low = T_low,P_nom_des=P_gro,T_amb_des=T_amb_des, N_cool=N_cooler) annotation(
      Placement(visible = true, transformation(origin = {-78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Turbine turbine(PR = PR, N_shaft = N_shaft, eta_design=eta_turb) annotation(
      Placement(visible = true, transformation(origin = {66, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.Exchanger exchanger(redeclare package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph,P_nom_des=P_gro,T_out_CO2_des=T_high,N_exch=N_exch, ratio_m_des=1) annotation(
      Placement(visible = true, transformation(origin = {46, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.CompressorOnShaft reCompressor(N_design = N_shaft,P_nom_des=P_gro,p_high_des=p_high) annotation(
      Placement(visible = true, transformation(origin = {-54, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.HeatRecuperatorDTAve LTR(N_q = N_LTR,P_nom_des=P_gro, ratio_m_des=1-gamma) annotation(
      Placement(visible = true, transformation(origin = {-42, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.FlowMixer mixer annotation(
      Placement(visible = true, transformation(origin = {-20, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.FlowSplitter splitter(gamma = gamma) annotation(
      Placement(visible = true, transformation(origin = {-62, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    
    parameter MedRec.ThermodynamicState state_HTF_in_des=MedRec.setState_pTX(1.0325*10^5,T_HTF_in_des);
    
     SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = 800+273.15, p_out = 10 ^ 5, m_flow = exchanger.m_HTF_des, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
        Placement(visible = true, transformation(origin = {-52, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {-50, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
    //Modelica.Blocks.Interfaces.RealInput parasities_internal;
    protected
    
  initial equation
   
    exchanger.state_HTF_des[N_exch] = MedRec.setState_pTX(10^5,T_high);
    exchanger.m_HTF_des=m_HTF_des;
    P_nom = (-turbine.W_turb_des) - mainCompressor.W_comp_des -reCompressor.W_comp_des- cooler.P_cool_des;
  
  // Thermodynamic state equalities
    //main loop
    exchanger.state_CO2_des[1] = HTR.state_comp_des[N_HTR];
    turbine.state_in_des = exchanger.state_CO2_des[N_exch];
    HTR.state_turb_des[N_HTR] = turbine.state_out_des;
    LTR.state_turb_des[N_LTR]=HTR.state_turb_des[1];
    cooler.state_des[1] = LTR.state_turb_des[1];
    mainCompressor.state_in_des = cooler.state_des[N_cooler];
    LTR.state_comp_des[1] = mainCompressor.state_out_des;
    // recompression loop
    reCompressor.state_in_des=LTR.state_turb_des[1];
    HTR.state_comp_des[1]=MedPB.setState_phX(p_high,ratio_m_des*LTR.state_comp_des[N_LTR].h+(1-ratio_m_des)*reCompressor.state_out_des.h);
  
  //mass flow equalities
    //main loop
    //exchanger.m_CO2_des = HTR.m_comp_des;
    turbine.m_des = exchanger.m_CO2_des;
    HTR.m_turb_des = turbine.m_des;
    LTR.m_turb_des=HTR.m_turb_des;
    cooler.m_des = LTR.m_turb_des*ratio_m_des;
    mainCompressor.m_des = cooler.m_des;
    LTR.m_comp_des=mainCompressor.m_des;
    //recompression loop
    HTR.m_comp_des = reCompressor.m_des+LTR.m_comp_des;
    reCompressor.m_des=gamma*LTR.m_turb_des;
    
    // Financial Analysis
    C_HTR= if MedPB.temperature(turbine.state_out_des)>=550+273.15 then 49.45*HTR.UA_HTR^0.7544*(1+0.02141*(MedPB.temperature(turbine.state_out_des)-550-273.15)) else 49.45*HTR.UA_HTR^0.7544;
    C_LTR=49.45*LTR.UA_HTR^0.7544;
    C_turbine= if TIT_CO2>= 550+273.15 then 406200*(-turbine.W_turb_des/10^6)^0.8*(1+1.137*10^(-5)*(TIT_CO2-550-273.15)^2) else 406200*(-turbine.W_turb_des/10^6)^0.8;
    C_mainCompressor = 1230000*(mainCompressor.W_comp_des/10^6)^0.3392;
    C_reCompressor = 1230000*(reCompressor.W_comp_des/10^6)^0.3392;
    C_cooler = 32.88*cooler.UA_cooler^0.75;
    C_generator = 108900*(P_nom/10^6)^0.5463;
    C_exchanger = pri_exchanger*exchanger.Q_HX_des*m_HTF_des/1000;
    C_PB=(C_HTR+C_LTR+C_turbine+C_mainCompressor+C_reCompressor+C_generator+C_cooler+C_exchanger)*1.05;
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
    //if external_parasities then
  //    connect(parasities_internal,parasities);
  //  else
  //    parasities_internal=0;
  //  end if;
    connect(LTR.from_turb_port_b, splitter.port_a) annotation(
      Line(points = {{-50, -34}, {-54, -34}, {-54, -36}, {-54, -36}}, color = {0, 127, 255}));
    connect(LTR.from_turb_port_a, HTR.from_turb_port_b) annotation(
      Line(points = {{-36, -34}, {0, -34}, {0, -32}, {0, -32}}, color = {0, 127, 255}));
    connect(reCompressor.port_b, mixer.second_port_a) annotation(
      Line(points = {{-50, -4}, {-20, -4}, {-20, -10}, {-20, -10}}, color = {0, 127, 255}));
    connect(LTR.from_comp_port_b, mixer.first_port_a) annotation(
      Line(points = {{-36, -24}, {-28, -24}, {-28, -18}, {-28, -18}}, color = {0, 127, 255}));
    connect(mixer.port_b, HTR.from_comp_port_a) annotation(
      Line(points = {{-12, -18}, {-4, -18}, {-4, -16}, {0, -16}}, color = {0, 127, 255}));
    connect(splitter.gamma_port_b, reCompressor.port_a) annotation(
      Line(points = {{-62, -28}, {-62, -28}, {-62, 8}, {-62, 8}}, color = {0, 127, 255}));
    connect(mainCompressor.port_b, LTR.from_comp_port_a) annotation(
      Line(points = {{-70, -18}, {-50, -18}, {-50, -24}, {-50, -24}}, color = {0, 127, 255}));
    connect(splitter.one_gamma_port_b, cooler.port_a) annotation(
      Line(points = {{-70, -36}, {-70, -36}, {-70, -62}, {-78, -62}, {-78, -62}}, color = {0, 127, 255}));
    connect(exchanger.CO2_port_a, HTR.from_comp_port_b) annotation(
      Line(points = {{38.8, 0.4}, {30.8, 0.4}, {30.8, 0.4}, {22.8, 0.4}, {22.8, -13.6}, {21.8, -13.6}, {21.8, -13.6}, {22.8, -13.6}}, color = {0, 127, 255}));
    connect(exchanger.CO2_port_b, turbine.port_a) annotation(
      Line(points = {{52.2, 0.2}, {59.2, 0.2}, {59.2, 0.2}, {64.2, 0.2}, {64.2, -23.8}, {62.2, -23.8}, {62.2, -24.8}, {62.2, -24.8}, {62.2, -25.8}}, color = {0, 127, 255}));
    connect(turbine.port_b, HTR.from_turb_port_a) annotation(
      Line(points = {{72, -36.6}, {48, -36.6}, {48, -34.6}, {22, -34.6}, {22, -34.1}, {22, -34.1}, {22, -31.6}}, color = {0, 127, 255}));
    connect(cooler.port_b, mainCompressor.port_a) annotation(
      Line(points = {{-78, -46}, {-88, -46}, {-88, -6}, {-82, -6}, {-82, -6}}, color = {0, 127, 255}));
    connect(cooler.T_amb, T_amb);
    connect(m_sup, exchanger.m_sup);
    
  
    m_sup = exchanger.HTF_port_a.m_flow >= exchanger.m_HTF_des * nu_min;
  //Closure equation
  if m_sup then 
    turbine.p_out=mainCompressor.p_out/PR;
  else
    exchanger.CO2_port_a.m_flow=exchanger.m_CO2_des;
  end if;
    
    eta_cycle = W_net / exchanger.Q_HX;
    der(E_net)=W_net;
    W_net = if m_sup then (-turbine.W_turb) - mainCompressor.W_comp - reCompressor.W_comp - cooler.P_cooling else 0;
    
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    annotation(
      Documentation(info = "<html>
  		<p>Direct-design model of a recompression sCO2 cycle with exchanger. On-design is performed in the initial equation section and off-design in the equation section.</p>
  <p>The mass flow in off-design is imposed by the turbine. The turbomachinery has the same set of variables: N, m_flow, p_out, h_out, with closure equation for 2 (performance maps). </p>
  <p> Therefore, one has to provide values for the two others. 
  <ul>
  <li>Compressors: N, m_flow provided. Outlet: p_out, h_out.</li>
  <li>Turbine: N, p_out provided. Outlet: m_flow, h_out. </li></ul>
   </p>
   <p> p_out,turb=p_out,comp/PR : this closure equation appeared to be the most logical, as one can interfere with a power block at full-load with the outlet pressure of the turbine, and that keeping the PR constant often augments the efficiency. It obviously has to be discussed. </p>
  
  <p> The currency is 2019$. The price of the exchanger is taken at 150$/kW_th because it is defined as the objective to reach for next-Gen CSP with particles</p>
  
  <p>N. T. Weiland, B. W. Lance, et S. R. Pidaparti, « SCO2 power cycle components cost correlations from DOE data spanning multiple scales and application », p. 17.</p>
  <p> Available at https://www.netl.doe.gov/projects/files/sCO2PowerCycleComponentCostCorrelationsfromDOEDataSpanningMultipleScalesandApplications_061819.pdf </p>
  <p>&nbsp;</p>
  		</html>"));
  end recompPB;

  model FlowSplitter
    extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedRec = SolarTherm.Media.CO2.CO2_ph;
    parameter Real gamma;
    Real gamma_var;
    Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b gamma_port_b(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {0, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {2, 84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b one_gamma_port_b(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    gamma_var = gamma;
    gamma_port_b.m_flow = -gamma_var * port_a.m_flow;
    one_gamma_port_b.m_flow = -(1 - gamma_var) * port_a.m_flow;
    gamma_port_b.p = port_a.p;
    one_gamma_port_b.p = port_a.p;
    gamma_port_b.h_outflow = inStream(port_a.h_outflow);
    one_gamma_port_b.h_outflow = inStream(port_a.h_outflow);
    port_a.h_outflow = inStream(gamma_port_b.h_outflow);
    annotation(
      Icon(graphics = {Text(origin = {0, 10}, extent = {{-56, -16}, {56, 16}}, textString = "SPLITTER")}),
      Diagram(graphics = {Text(origin = {7, 8}, extent = {{-49, -16}, {49, 16}}, textString = "SPLITTER")}));
  end FlowSplitter;

  model FlowMixer "This model is useful for the recompression cycle cycle, as it allows to mix two different fluid. The pressure in both should entrance should be the same; in case it is not, we ponderated it by the mass flows: as it is the same molar mass, the resulting pressure should look like that."
    extends SolarTherm.Media.CO2.PropCO2;
    replaceable package MedRec = SolarTherm.Media.CO2.CO2_ph;
    Modelica.Fluid.Interfaces.FluidPort_a first_port_a(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {0, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_a second_port_a(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {2, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedRec) annotation(
      Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    port_b.m_flow = -(first_port_a.m_flow + second_port_a.m_flow);
    port_b.h_outflow = (first_port_a.m_flow * inStream(first_port_a.h_outflow) + second_port_a.m_flow * inStream(second_port_a.h_outflow)) / (first_port_a.m_flow + second_port_a.m_flow);
    port_b.p = (first_port_a.m_flow * first_port_a.p + second_port_a.m_flow * second_port_a.p) / (first_port_a.m_flow + second_port_a.m_flow);
    first_port_a.h_outflow = inStream(port_b.h_outflow);
    second_port_a.h_outflow = inStream(port_b.h_outflow);
    annotation(
      Diagram(graphics = {Text(origin = {3, 13}, extent = {{-39, -17}, {39, 17}}, textString = "MIXER")}),
      Icon(graphics = {Text(origin = {6, 14}, extent = {{-44, -28}, {44, 28}}, textString = "MIXER")}));
  end FlowMixer;

  model testHX
    extends SolarTherm.Media.CO2.PropCO2;
    parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
    parameter Real m_des = 100;
    parameter Real m_flow = 90;
    SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow srcTLMDt(p_out = 230 * 10 ^ 5, T_out = 550, m_flow = m_flow) annotation(
      Placement(visible = true, transformation(origin = {46, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    SourceFlow srcTLMDc(p_out = 85 * 10 ^ 5, T_out = 120, m_flow = 0.7 * m_flow) annotation(
      Placement(visible = true, transformation(origin = {-44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sinkTLMDt annotation(
      Placement(visible = true, transformation(origin = {-44, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sinkTLMDc annotation(
      Placement(visible = true, transformation(origin = {42, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    SourceFlow srcDTAvec(p_out = 85 * 10 ^ 5, T_out = 120, m_flow = 0.7 * m_flow) annotation(
      Placement(visible = true, transformation(origin = {-42, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow srcDTAvet(p_out = 230 * 10 ^ 5, T_out = 550, m_flow = m_flow) annotation(
      Placement(visible = true, transformation(origin = {46, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    SinkFlow sinkDTAvet annotation(
      Placement(visible = true, transformation(origin = {-42, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sinkDTAvec annotation(
      Placement(visible = true, transformation(origin = {46, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.DirectDesign.HeatRecuperatorDTAve TLMD(N_q = 15, P_nom_des = 10 ^ 7) annotation(
      Placement(visible = true, transformation(origin = {0, -16}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    DirectDesign.HeatRecuperatorDTAve DTAve(N_q = 60, P_nom_des = 10 ^ 7) annotation(
      Placement(visible = true, transformation(origin = {-4, 48}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  initial equation
    DTAve.h_in_comp_des = 500000;
    DTAve.h_in_turb_des = 10 ^ 6;
    DTAve.p_in_comp_des = 250 * 10 ^ 5;
    DTAve.p_in_turb_des = 80 * 10 ^ 5;
    DTAve.m_comp_des = 0.7 * m_des;
    DTAve.m_turb_des = m_des;
    TLMD.h_in_comp_des = 500000;
    TLMD.h_in_turb_des = 10 ^ 6;
    TLMD.p_in_comp_des = 250 * 10 ^ 5;
    TLMD.p_in_turb_des = 85 * 10 ^ 5;
    TLMD.m_comp_des = 0.7 * m_des;
    TLMD.m_turb_des = m_des;
  equation
    connect(DTAve.from_comp_port_b, sinkDTAvec.port_a) annotation(
      Line(points = {{6, 56}, {6, 56}, {6, 66}, {38, 66}, {38, 66}}, color = {0, 127, 255}));
    connect(srcDTAvec.port_b, DTAve.from_comp_port_a) annotation(
      Line(points = {{-34, 64}, {-16, 64}, {-16, 54}, {-16, 54}}, color = {0, 127, 255}));
    connect(DTAve.from_turb_port_b, sinkDTAvet.port_a) annotation(
      Line(points = {{-16, 38}, {-14, 38}, {-14, 30}, {-34, 30}, {-34, 30}}, color = {0, 127, 255}));
    connect(DTAve.from_turb_port_a, srcDTAvet.port_b) annotation(
      Line(points = {{6, 38}, {6, 38}, {6, 30}, {38, 30}, {38, 30}}, color = {0, 127, 255}));
    connect(TLMD.from_comp_port_a, srcTLMDc.port_b) annotation(
      Line(points = {{-14, -8}, {-16, -8}, {-16, 0}, {-36, 0}, {-36, 0}}, color = {0, 127, 255}));
    connect(TLMD.from_turb_port_b, sinkTLMDt.port_a) annotation(
      Line(points = {{-12, -26}, {-14, -26}, {-14, -30}, {-36, -30}, {-36, -30}}, color = {0, 127, 255}));
    connect(TLMD.from_turb_port_a, srcTLMDt.port_b) annotation(
      Line(points = {{12, -26}, {38, -26}, {38, -30}, {38, -30}}, color = {0, 127, 255}));
    connect(TLMD.from_comp_port_b, sinkTLMDc.port_a) annotation(
      Line(points = {{12, -8}, {12, -8}, {12, 0}, {34, 0}, {34, 0}}, color = {0, 127, 255}));
  end testHX;
  annotation(
      Documentation(info = "<html>
  		<p>This section proposes direct-design calculation of sCO2 cycles. The on-design is performed in the initial equation section, and off-design is performed in the equation section. </p>
  		<p> It is therefore necessary to connect cycles between them in the initial section as well, through thermodynamic states at the inlet and outlet. </p>
  		<p> It was developped by the CEA (France) from the thesis of J.J. Dyreby, MIT. It still needs validation. </p> 
  		<p> Calculation times are a bit heavy; this should be used for validation of the optimization performed with QuickCO2PB </p>
  		</html>"));
end DirectDesign;