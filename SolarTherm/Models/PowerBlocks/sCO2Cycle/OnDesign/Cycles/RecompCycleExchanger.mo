within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model RecompCycleExchanger "On-design calculation of the sCO2 recompression cycle. Check the code of the SourceFlow for errors."
      import FI = SolarTherm.Models.Analysis.Finances;
      import SI = Modelica.SIunits;
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      
      inner Modelica.Fluid.System system(
    		energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial,
    allowFlowReversal=false);
    
      
      //Parameters
      parameter SI.AbsolutePressure p_high = 250 * 10 ^ 5 "high pressure of the cycle";
      parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_HTF_in = 800 + 273.15 "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_amb = 30 + 273.15 "ambiant  temperature";
      parameter Real PR = 2.777 "Pressure ratio";
      parameter SI.Power P_nom =10*10^6 "Nominal electrical power";
      
      parameter SI.Efficiency eta_comp = 0.89 "Isentropic efficiency of the compressors";
      parameter SI.Efficiency eta_turb = 0.93 "Isentropic efficiency of the turbine";
      parameter Real gamma = 0.284 "Part of the mass flow going to the recompression directly";
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
      parameter Integer N_exch = 5 "Number of discretization of the exchanger";
      parameter Real pinchRecuperator = 5;
      parameter Real[4] choiceN = {75000,30000,10000,3600}*5/6*0.10471975512 ;
      parameter Real N_shaft=(choiceN[integer(Modelica.Math.log(P_nom/10^6)/Modelica.Math.log(10))+2]);
      
      // Financial prices
      
     
      parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
       
      // Instanciation of the components
      
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine turbine(PR = PR, eta_turb = eta_turb, T_amb = T_amb,N_shaft=N_shaft) annotation(
        Placement(visible = true, transformation(origin = {74, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft mainCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb,N_compressor=N_shaft) annotation(
        Placement(visible = true, transformation(origin = {-60, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler(T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-66, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft reCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb,N_compressor=N_shaft) annotation(
        Placement(visible = true, transformation(origin = {-24, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve LTR(N_q = N_q, flowGuess = 700, T_amb = T_amb, pinchRecuperator = pinchRecuperator) annotation(
        Placement(visible = true, transformation(origin = {-26, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTR(N_q = N_q, flowGuess = 700, T_amb = T_amb, pinchRecuperator = pinchRecuperator) annotation(
        Placement(visible = true, transformation(origin = {26, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowSplitter splitter(gamma = gamma) annotation(
        Placement(visible = true, transformation(origin = {-50, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowMixer mixer annotation(
        Placement(visible = true, transformation(origin = {-2, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Exchanger exchanger (T_amb=T_amb, TIT_CO2=T_high, redeclare package MedRec=SolarTherm.Media.SolidParticles.CarboHSP_ph) annotation(
        Placement(visible = true, transformation(origin = {46, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
      // Variables to investigate the cycle and its simulation.
      
      SI.Efficiency efficiencyCycle "Efficiency of the cycle";
      Real E_bal_check;
      SolarTherm.Types.SpecificWork W_out "Specific Work of the cycle";
      
      //Financial Analysis
      FI.Money C_LTR "cost of the low-temperature recuperator";
      FI.Money C_HTR "cost of the high-temperature recuperator";
      FI.Money C_turbine "cost of the turbine";
      FI.Money C_mainCompressor "cost of the compressor";
      FI.Money C_reCompressor "cost of the reCompressor";
      FI.Money C_exchanger "cost of the exchanger";
      FI.Money C_generator "cost of the generator";
      FI.Money C_cooler "cost of the cooler";
      FI.Money C_PB "Overall cost of the power block";
    
      SI.Efficiency eta_carnot;
      // Exergy analysis
      Real ex_d_percent_mainCompressor "MainCompressor exergy destruction";
      Real ex_d_percent_LTR "LTR exergy destruction";
      Real ex_d_percent_HTRecuperator "HTRecuperator exergy destruction";
      Real ex_d_percent_reCompressor "reCompressor exergy destruction";
      Real ex_d_percent_exchanger "heater exergy destruction";
      Real ex_d_percent_turbine "turbine exergy destruction";
      Real ex_d_percent_cooler "cooler exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_tot "Total exergy destruction";
      SolarTherm.Types.SpecificExergy ex_in "Inlet of exergy at the heater";
      SI.Efficiency eta_ex "Exergetic efficiency = P_nom/ex_in";
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = T_HTF_in, use_m_parameter = false, m_flow =730, p_out = 10 ^ 5, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph) annotation(
        Placement(visible = true, transformation(origin = {74, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {22, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(src.port_b, exchanger.HTF_port_a) annotation(
        Line(points = {{66, 30}, {52, 30}, {52, 26}, {52, 26}}, color = {0, 127, 255}));
      connect(sink.port_a, exchanger.HTF_port_b) annotation(
        Line(points = {{30, 28}, {38, 28}, {38, 26}, {38, 26}}, color = {0, 127, 255}));
      connect(exchanger.CO2_port_b, turbine.port_a) annotation(
        Line(points = {{52, 16}, {70, 16}, {70, 14}, {70, 14}}, color = {0, 127, 255}));
      connect(HTR.from_comp_port_b, exchanger.CO2_port_a) annotation(
        Line(points = {{32, -8}, {38, -8}, {38, 16}, {38, 16}}, color = {0, 127, 255}));
// Connectors
      connect(splitter.one_gamma_port_b, cooler.port_a) annotation(
        Line(points = {{-58, -18}, {-58, -18}, {-58, -22}, {-66, -22}, {-66, -22}}, color = {0, 127, 255}));
      connect(splitter.port_a, LTR.from_turb_port_b) annotation(
        Line(points = {{-42, -18}, {-34, -18}, {-34, -18}, {-34, -18}}, color = {0, 127, 255}));
      connect(splitter.gamma_port_b, reCompressor.port_a) annotation(
        Line(points = {{-50, -10}, {-48, -10}, {-48, 20}, {-32, 20}, {-32, 20}}, color = {0, 127, 255}));
      connect(mainCompressor.port_b, LTR.from_comp_port_a) annotation(
        Line(points = {{-56, 10}, {-34, 10}, {-34, -8}, {-34, -8}, {-34, -8}}, color = {0, 127, 255}));
      connect(HTR.from_turb_port_b, LTR.from_turb_port_a) annotation(
        Line(points = {{18, -18}, {-20, -18}, {-20, -18}, {-20, -18}}, color = {0, 127, 255}));
      connect(LTR.from_comp_port_b, mixer.first_port_a) annotation(
        Line(points = {{-20, -8}, {-10, -8}, {-10, -8}, {-10, -8}}, color = {0, 127, 255}));
      connect(reCompressor.port_b, mixer.second_port_a) annotation(
        Line(points = {{-20, 8}, {0, 8}, {0, 0}, {-2, 0}, {-2, 0}}, color = {0, 127, 255}));
      connect(HTR.from_comp_port_a, mixer.port_b) annotation(
        Line(points = {{18, -8}, {6, -8}, {6, -8}, {6, -8}}, color = {0, 127, 255}));
      
      connect(turbine.port_b, HTR.from_turb_port_a) annotation(
        Line(points = {{80, 4}, {80, 4}, {80, -18}, {32, -18}, {32, -18}}, color = {0, 127, 255}));
      connect(cooler.port_b, mainCompressor.port_a) annotation(
        Line(points = {{-66, -6}, {-80, -6}, {-80, 22}, {-68, 22}, {-68, 22}, {-68, 22}}, color = {0, 127, 255}));
     

// Cycle efficiency. Mass flow is imposed by the power equation.
      P_nom = (-turbine.W_turb) - (mainCompressor.W_comp + reCompressor.W_comp)-cooler.P_cooling;
    
      efficiencyCycle * exchanger.Q_HX = P_nom;
      E_bal_check = turbine.W_turb + mainCompressor.W_comp + reCompressor.W_comp + exchanger.Q_HX + cooler.Q_cooler;
      W_out = P_nom / turbine.port_a.m_flow;

    // Financial analysis
    
      C_LTR = 49.45*LTR.UA_HTR^0.7544;
      C_HTR= if HTR.T_from_turb[N_q]>=550+273.15 then 49.45*HTR.UA_HTR^0.7544*(1+0.02141*(HTR.T_from_turb[N_q]-550-273.15)) else 49.45*HTR.UA_HTR^0.7544;
      C_turbine= if T_high>= 550+273.15 then 406200*(-turbine.W_turb/10^6)^0.8*(1+1.137*10^(-5)*(T_high-550-273.15)^2) else 406200*(-turbine.W_turb/10^6)^0.8;
      C_mainCompressor = 1230000*(mainCompressor.W_comp/10^6)^0.3392;
      C_reCompressor = 1230000*(reCompressor.W_comp/10^6)^0.3392;
      C_cooler = 32.88*cooler.UA_cooler^0.75;
      C_generator = 108900*(P_nom/10^6)^0.5463;
      C_exchanger = pri_exchanger*exchanger.Q_HX/1000;
      C_PB=(C_HTR+C_LTR+C_turbine+C_mainCompressor+C_reCompressor+C_generator+C_cooler+C_exchanger)*1.05;
      // 1.05 corresponds to inflation from 2017 to 2019, as correlations are in 2017' U.S. dollars.
  // Exergy analysis
    
      ex_d_tot = mainCompressor.ex_d + LTR.ex_d + HTR.ex_d + reCompressor.ex_d + exchanger.ex_d + turbine.ex_d + cooler.ex_d;
      mainCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_mainCompressor;
      LTR.ex_d * 100 / ex_d_tot = ex_d_percent_LTR;
      HTR.ex_d * 100 / ex_d_tot = ex_d_percent_HTRecuperator;
      reCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_reCompressor;
      exchanger.ex_d * 100 / ex_d_tot = ex_d_percent_exchanger;
      turbine.ex_d * 100 / ex_d_tot = ex_d_percent_turbine;
      cooler.ex_d * 100 / ex_d_tot = ex_d_percent_cooler;
      ex_in = exchanger.Q_HX * (1 - T_amb / T_high);
      eta_ex = 1 - ex_d_tot / ex_in;
      eta_carnot = 1 - T_amb / T_high;
      annotation(
      Icon,
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
  __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
      annotation(
        Documentation(info = "<html>
    		<p>On-design model of a recompression sCO2 cycle with exchanger. A number of discretization of the heat recuperators of 15 seems accurate when comparated with 60, regarding numerical complexity.</p>
    <p>The mass flow is either imposed by the nominal power or by the HTF mass flow rate. An equation on the mass flow rate of the CO2 and the HTF is necessary; it was taken such as m_flow*Cp would be even for both flows.</p>
    <p>A calculation of the price is performed based on a cost estimation of the different components, from  Weiland et al.The uncertainty is between -30%/50%. Depending on the power block layout chosen, other correlations might have to be implemented (motor for the compressor, gearbox, ..). See the article for more informations.</p>
    <p> The currency is 2017$. The price of the exchanger is taken at 150$/kW_th because it is defined as the objective to reach for next-Gen CSP with particles</p>
    <p>An exergy analysis is implemented based on a class from Pr. Neveu (UPVD).</p>
    <p>N. T. Weiland, B. W. Lance, et S. R. Pidaparti, « SCO2 power cycle components cost correlations from DOE data spanning multiple scales and application », p. 17.</p>
    <p> Available at https://www.netl.doe.gov/projects/files/sCO2PowerCycleComponentCostCorrelationsfromDOEDataSpanningMultipleScalesandApplications_061819.pdf </p>
    <p>&nbsp;</p>
    		</html>"));
    end RecompCycleExchanger;