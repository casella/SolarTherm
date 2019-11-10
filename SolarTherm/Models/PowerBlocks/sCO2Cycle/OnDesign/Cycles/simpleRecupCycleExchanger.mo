within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model simpleRecupCycleExchanger "Simple recuperation cycle with an exchanger"
      extends SolarTherm.Media.CO2.PropCO2;
      import FI = FI;
      import SI = Modelica.SIunits;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      //Parameters
      parameter SI.AbsolutePressure p_high = 250 * 10 ^ 5 "high pressure of the cycle";
      parameter SI.ThermodynamicTemperature T_high = 800 + 273.15 "inlet temperature of the HTF";
      parameter SI.ThermodynamicTemperature T_amb = 30 + 273.15 "ambiant temperature";
      parameter Real PR = 2.5 "Pressure ratio";
      parameter SI.Power P_nom = 100*10^6 "Nominal electrical power";
      parameter SI.Efficiency eta_comp = 0.89 "Isentropic efficiency of the compressors";
      parameter SI.Efficiency eta_turb = 0.93 "Isentropic efficiency of the turbine";
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
    
      // Instanciation of the components
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine turbine(PR = PR, eta_turb = eta_turb, T_amb = T_amb, is_second_turbine = true) annotation(
        Placement(visible = true, transformation(origin = {74, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft mainCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-60, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler(T_amb = T_amb, P_nom = P_nom) annotation(
        Placement(visible = true, transformation(origin = {-66, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTR(N_q = N_q, flowGuess = 120, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-2, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Exchanger exchanger(pinch = 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {36, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //equivalent to T_Sodium = 700K
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sinkMS annotation(
        Placement(visible = true, transformation(origin = {-2, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Variables to investigate the cycle and its simulation.
      SI.Efficiency efficiencyCycle;
      Real E_bal_check;
      FI.Money C_recuperators "cost of the recuperators";
      FI.Money C_turbomachinery "cost of the turbomachinery";
      FI.Money C_exchanger "cost of the exchanger";
      FI.Money C_cooler "cost of the cooler";
      FI.Money C_PB "Overall cost of the power block";
      //Exergy analysis
      Real ex_d_percent_Compressor "Compressor exergy destruction";
      Real ex_d_percent_HTR "Recuperator exergy destruction";
      Real ex_d_percent_exchanger "heater exergy destruction";
      Real ex_d_percent_turbine "turbine exergy destruction";
      Real ex_d_percent_cooler "cooler exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_tot "Total exergy destruction";
      SolarTherm.Types.SpecificExergy ex_in;
      SI.Efficiency eta_ex "Exergetic efficiency";
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow srcMS(T_out = T_high - 273.15, use_m_parameter = false, m_flow = 1.48, p_out = 10 ^ 5, redeclare package MedPB = SolarTherm.Media.Sodium.ConstSodium) annotation(
        Placement(visible = true, transformation(origin = {60, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    equation
      connect(srcMS.port_b, exchanger.HTF_port_a) annotation(
        Line(points = {{52, 40}, {40, 40}, {40, 16}, {42, 16}, {42, 16}}, color = {0, 127, 255}));
      connect(sinkMS.port_a, exchanger.HTF_port_b) annotation(
        Line(points = {{6, 30}, {28, 30}, {28, 16}, {28, 16}, {28, 16}}, color = {0, 127, 255}));
//Connections
      connect(exchanger.CO2_port_b, turbine.port_a) annotation(
        Line(points = {{42, 6}, {58, 6}, {58, 14}, {70, 14}, {70, 14}}, color = {0, 127, 255}));
      connect(HTR.from_comp_port_b, exchanger.CO2_port_a) annotation(
        Line(points = {{4, -4}, {28, -4}, {28, 6}, {28, 6}, {28, 6}}, color = {0, 127, 255}));
      connect(mainCompressor.port_b, HTR.from_comp_port_a) annotation(
        Line(points = {{-56, 10}, {-10, 10}, {-10, -4}}, color = {0, 127, 255}));
      connect(HTR.from_turb_port_b, cooler.port_a) annotation(
        Line(points = {{-9, -14}, {-9, -22}, {-66, -22}}, color = {0, 127, 255}));
      connect(HTR.from_turb_port_a, turbine.port_b) annotation(
        Line(points = {{4, -14}, {80, -14}, {80, 4}}, color = {0, 127, 255}));
      connect(cooler.port_b, mainCompressor.port_a) annotation(
        Line(points = {{-66, -6}, {-80, -6}, {-80, 22}, {-68, 22}, {-68, 22}, {-68, 22}}, color = {0, 127, 255}));
//fixes the mass flow
  P_nom = (-turbine.W_turb) - mainCompressor.W_comp-cooler.P_cooling;
//Calculates the energetic efficiency of the cycle
      efficiencyCycle * exchanger.Q_HX_HTR = P_nom;
//check the 1st law of thermodynamics
      E_bal_check = turbine.W_turb + mainCompressor.W_comp + exchanger.Q_HX_HTR + cooler.Q_cooler;
      
      //Financial analysis
      C_recuperators = pri_recuperators*(HTR.UA_HTR)/1000;
      C_turbomachinery = pri_turbomachinery*P_nom/1000;
      C_cooler = pri_cooler*cooler.UA_cooler/1000;
      C_exchanger = pri_exchanger*exchanger.UA_TLMD/1000;
      C_PB=C_recuperators+C_turbomachinery+C_cooler+C_exchanger;
      
// Exergy analysis
      ex_d_tot = mainCompressor.ex_d + HTR.ex_d + exchanger.ex_d + turbine.ex_d + cooler.ex_d;
      mainCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_Compressor;
      HTR.ex_d * 100 / ex_d_tot = ex_d_percent_HTR;
      exchanger.ex_d * 100 / ex_d_tot = ex_d_percent_exchanger;
      turbine.ex_d * 100 / ex_d_tot = ex_d_percent_turbine;
      cooler.ex_d * 100 / ex_d_tot = ex_d_percent_cooler;
      ex_in = exchanger.Q_HX_HTR * (1 - T_amb / T_high);
      eta_ex = ((-turbine.W_turb) - mainCompressor.W_comp) / ex_in;
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
        Documentation(info = "<html>
    		<p>On-design model of a simple recuperation sCO2 cycle with an HTF/sCO2 counter-flow HX. A number of discretization of the heat recuperators of 15 seems accurate when compared with 60, regarding numerical complexity.</p>
    
    <p>A calculation of the price is performed based on a cost estimation of the different components, from  Weiland et al.The uncertainty is between -30%/50%. Depending on the power block layout chosen, other correlations might have to be implemented (motor for the compressor, gearbox, ..). See the article for more informations.</p>
    <p> The currency is 2017$, except for the overall PB where it is expressed in 2019 US$. The price of the heater is taken at 150$/kW_th because it is defined as the objective to reach for next-Gen CSP with particles</p>
    <p>An exergy analysis is implemented based on a class from Pr. Neveu (UPVD).</p>
    <p>N. T. Weiland, B. W. Lance, et S. R. Pidaparti, « SCO2 power cycle components cost correlations from DOE data spanning multiple scales and application », p. 17.</p>
    <p> Available at https://www.netl.doe.gov/projects/files/sCO2PowerCycleComponentCostCorrelationsfromDOEDataSpanningMultipleScalesandApplications_061819.pdf </p>
    <p>&nbsp;</p>
    		</html>"));
    end simpleRecupCycleExchanger;