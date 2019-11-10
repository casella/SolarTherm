within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model simpleRecupCycleHeater
      import FI = SolarTherm.Models.Analysis.Finances;
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      //Parameters
      parameter SI.AbsolutePressure p_high = CV.from_bar(200) "high pressure of the cycle";
      parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_amb = CV.from_degC(40) "ambiant temperature";
      parameter SI.TemperatureDifference deltaT_cooler = 15 "Approach difference of temperature at the outlet of the cooler";
      parameter Real PR = 2.313 "Pressure ratio";
      parameter SI.Power P_nom = 100*10^6 "Nominal electrical power";
      parameter SI.Efficiency eta_comp = 0.87 "Isentropic efficiency of the compressors";
      parameter SI.Efficiency eta_turb = 0.9 "Isentropic efficiency of the turbine";
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
      // Instanciation of the components
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine turbine(PR = PR, eta_turb = eta_turb, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {74, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft mainCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-60, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler(T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-66, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTR(N_q = N_q, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-2, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Heater heater(T_high = T_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {30, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Variables to investigate the cycle and its simulation.
      SI.Efficiency efficiencyCycle;
      Real E_bal_check;
      //Exergy analysis
      SolarTherm.Types.SpecificExergy ex_d_percent_Compressor "Compressor exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_percent_HTR "Recuperator exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_percent_heater "heater exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_percent_turbine "turbine exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_percent_cooler "cooler exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_tot "Total exergy destruction";
      SolarTherm.Types.SpecificExergy ex_in;
      SI.Efficiency eta_ex "Exergetic efficiency";
      //Financial anaylisis
      FI.Money C_HTR  "cost of the heat recuperator";
      FI.Money C_turbine "cost of the turbine";
      FI.Money C_compressor "cost of the compressor";
      FI.Money C_exchanger "cost of the exchanger";
      FI.Money C_generator "cost of the generator";
      FI.Money C_cooler "cost of the cooler";
      FI.Money C_PB "Overall cost of the power block";
      FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
    equation
//Connections
      connect(heater.port_b, turbine.port_a) annotation(
        Line(points = {{30, 14}, {70, 14}, {70, 14}, {70, 14}}, color = {0, 127, 255}));
      connect(HTR.from_comp_port_b, heater.port_a) annotation(
        Line(points = {{4, -4}, {30, -4}, {30, -2}, {30, -2}}, color = {0, 127, 255}));
      connect(mainCompressor.port_b, HTR.from_comp_port_a) annotation(
        Line(points = {{-56, 10}, {-10, 10}, {-10, -4}}, color = {0, 127, 255}));
      connect(HTR.from_turb_port_b, cooler.port_a) annotation(
        Line(points = {{-9, -14}, {-9, -22}, {-66, -22}}, color = {0, 127, 255}));
      connect(HTR.from_turb_port_a, turbine.port_b) annotation(
        Line(points = {{4, -14}, {80, -14}, {80, 4}}, color = {0, 127, 255}));
      connect(cooler.port_b, mainCompressor.port_a) annotation(
        Line(points = {{-66, -6}, {-80, -6}, {-80, 22}, {-68, 22}, {-68, 22}, {-68, 22}}, color = {0, 127, 255}));
//fixes the mass flow
      P_nom = (-turbine.W_turb) - mainCompressor.W_comp;
//Calculates the energetic efficiency of the cycle
      efficiencyCycle * heater.Q_heater = P_nom;
//check the 1st law of thermodynamics
      E_bal_check = turbine.W_turb + mainCompressor.W_comp + heater.Q_heater + cooler.Q_cooler;
// Exergy analysis
      ex_d_tot = mainCompressor.ex_d + HTR.ex_d + heater.ex_d + turbine.ex_d + cooler.ex_d;
      mainCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_Compressor;
      HTR.ex_d * 100 / ex_d_tot = ex_d_percent_HTR;
      heater.ex_d * 100 / ex_d_tot = ex_d_percent_heater;
      turbine.ex_d * 100 / ex_d_tot = ex_d_percent_turbine;
      cooler.ex_d * 100 / ex_d_tot = ex_d_percent_cooler;
      ex_in = heater.Q_heater * (1 - T_amb / T_high);
      eta_ex = ((-turbine.W_turb) - mainCompressor.W_comp) / ex_in;
      
      //Financial analysis
      C_HTR= if HTR.T_from_turb[N_q]>=550+273.15 then 49.45*HTR.UA_HTR^0.7544*(1+0.02141*(HTR.T_from_turb[N_q]-550-273.15)) else 49.45*HTR.UA_HTR^0.7544;
      C_turbine= if T_high>= 550+273.15 then 406200*(-turbine.W_turb/10^6)^0.8*(1+1.137*10^(-5)*(T_high-550-273.15)^2) else 406200*(-turbine.W_turb/10^6)^0.8;
      C_compressor = 1230000*(mainCompressor.W_comp/10^6)^0.3392;
      C_cooler = 32.88*cooler.UA_cooler^0.75;
      C_generator = 108900*(P_nom/10^6)^0.5463;
      C_exchanger = pri_exchanger*heater.Q_heater/1000;
      C_PB=(C_HTR+C_turbine+C_compressor+C_generator+C_cooler+C_exchanger)*1.05;
      // 1.05 corresponds to inflation from 2017 to 2019, as correlations are in 2017' U.S. dollars.
      annotation(
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
        Documentation(info = "<html>
    		<p>On-design model of a simple recuperation sCO2 cycle with a simple heater. A number of discretization of the heat recuperators of 15 seems accurate when compared with 60, regarding numerical complexity.</p>
    
    <p>A calculation of the price is performed based on a cost estimation of the different components, from  Weiland et al.The uncertainty is between -30%/50%. Depending on the power block layout chosen, other correlations might have to be implemented (motor for the compressor, gearbox, ..). See the article for more informations.</p>
    <p> The currency is 2017$, except for the overall PB where it is expressed in 2019 US$. The price of the heater is taken at 150$/kW_th because it is defined as the objective to reach for next-Gen CSP with particles</p>
    <p>An exergy analysis is implemented based on a class from Pr. Neveu (UPVD).</p>
    <p>N. T. Weiland, B. W. Lance, et S. R. Pidaparti, « SCO2 power cycle components cost correlations from DOE data spanning multiple scales and application », p. 17.</p>
    <p> Available at https://www.netl.doe.gov/projects/files/sCO2PowerCycleComponentCostCorrelationsfromDOEDataSpanningMultipleScalesandApplications_061819.pdf </p>
    <p>&nbsp;</p>
    		</html>"));
    end simpleRecupCycleHeater;