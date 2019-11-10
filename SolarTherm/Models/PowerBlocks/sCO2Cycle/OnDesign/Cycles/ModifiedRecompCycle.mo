within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model ModifiedRecompCycle "Model of the modified sCO2 recompression cycle. "
      extends SolarTherm.Media.CO2.PropCO2;
      //Parameters
      parameter Modelica.SIunits.AbsolutePressure p_high = 250 * 10 ^ 5 "high pressure of the cycle";
      parameter Modelica.SIunits.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the turbine";
      parameter Modelica.SIunits.ThermodynamicTemperature T_amb = 40 + 273.15 "ambiant  temperature";
      parameter Real PR = 2.371 "Pressure ratio";
      parameter Modelica.SIunits.Power P_nom = 10 ^ 8 "Nominal electrical power";
      parameter Modelica.SIunits.Efficiency eta_comp = 0.87 "Isentropic efficiency of the compressors";
      parameter Modelica.SIunits.Efficiency eta_turb = 0.9 "Isentropic efficiency of the turbine";
      parameter Real gamma = 0.7 "Part of the mass flow going to the recompression directly";
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
      // Instanciation of the components
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Heater heater(T_high = T_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {32, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine turbine(PR = PR, eta_turb = eta_turb, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {74, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler(T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-66, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft reCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb, flowGuess = P_nom / 10 ^ 5 * (1 - gamma)) annotation(
        Placement(visible = true, transformation(origin = {2, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve LTRecuperator(N_q = N_q, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-26, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTRecuperator1(N_q = N_q, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {26, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowMixer mixer annotation(
        Placement(visible = true, transformation(origin = {-2, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Variables to investigate the cycle and its simulation.
      Modelica.SIunits.Efficiency efficiencyCycle "Efficiency of the cycle";
      Real E_bal_check;
      SolarTherm.Types.SpecificWork W_out "Specific Work of the cycle";
      SolarTherm.Types.Conductance UA_LTR(start = 2 * 10 ^ 6);
      SolarTherm.Types.Conductance UA_HTR(start = 2 * 10 ^ 6);
      //Modelica.SIunits.Efficiency eta_carnot;
      // Exergy analysis
      //  Real ex_d_percent_mainCompressor "MainCompressor exergy destruction";
      //  Real ex_d_percent_LTRecuperator "LTRecuperator exergy destruction";
      //  Real ex_d_percent_HTRecuperator "HTRecuperator exergy destruction";
      //  Real ex_d_percent_reCompressor "reCompressor exergy destruction";
      //  Real ex_d_percent_heater "heater exergy destruction";
      //  Real ex_d_percent_turbine "turbine exergy destruction";
      //  Real ex_d_percent_cooler "cooler exergy destruction";
      //  SolarTherm.Types.SpecificExergy ex_d_tot "Total exergy destruction";
      //  SolarTherm.Types.SpecificExergy ex_in "Inlet of exergy at the heater";
      //  Modelica.SIunits.Efficiency eta_ex "Exergetic efficiency = P_nom/ex_in";
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft firstCompressor(PR = PR, eta_comp = eta_comp, p_out = p_high, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-66, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowSplitter splitter(gamma = gamma) annotation(
        Placement(visible = true, transformation(origin = {-48, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler secondCooler annotation(
        Placement(visible = true, transformation(origin = {-34, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft secondCompressor annotation(
        Placement(visible = true, transformation(origin = {-34, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(secondCompressor.port_b, LTRecuperator.from_comp_port_a) annotation(
        Line(points = {{-30, 0}, {-34, 0}, {-34, -8}, {-34, -8}}, color = {0, 127, 255}));
      connect(secondCooler.port_b, secondCompressor.port_a) annotation(
        Line(points = {{-34, 20}, {-42, 20}, {-42, 12}, {-42, 12}, {-42, 12}}, color = {0, 127, 255}));
      connect(splitter.one_gamma_port_b, reCompressor.port_a) annotation(
        Line(points = {{-40, 48}, {-6, 48}, {-6, 24}, {-6, 24}, {-6, 24}}, color = {0, 127, 255}));
      connect(splitter.gamma_port_b, secondCooler.port_a) annotation(
        Line(points = {{-48, 40}, {-34, 40}, {-34, 36}, {-34, 36}, {-34, 36}}, color = {0, 127, 255}));
      connect(firstCompressor.port_b, splitter.port_a) annotation(
        Line(points = {{-58, 18}, {-58, 18}, {-58, 48}, {-56, 48}, {-56, 48}}, color = {0, 127, 255}));
      connect(cooler.port_b, firstCompressor.port_a) annotation(
        Line(points = {{-66, -6}, {-70, -6}, {-70, 6}, {-70, 6}, {-70, 6}}, color = {0, 127, 255}));
      connect(LTRecuperator.from_turb_port_b, cooler.port_a) annotation(
        Line(points = {{-34, -18}, {-34, -18}, {-34, -22}, {-66, -22}, {-66, -22}}, color = {0, 127, 255}));
      connect(reCompressor.port_b, mixer.second_port_a) annotation(
        Line(points = {{5, 13}, {0, 13}, {0, 0}, {-2, 0}}, color = {0, 127, 255}));
      connect(HTRecuperator1.from_turb_port_b, LTRecuperator.from_turb_port_a) annotation(
        Line(points = {{18, -18}, {-20, -18}, {-20, -18}, {-20, -18}}, color = {0, 127, 255}));
      connect(LTRecuperator.from_comp_port_b, mixer.first_port_a) annotation(
        Line(points = {{-20, -8}, {-10, -8}, {-10, -8}, {-10, -8}}, color = {0, 127, 255}));
      connect(HTRecuperator1.from_comp_port_a, mixer.port_b) annotation(
        Line(points = {{18, -8}, {6, -8}, {6, -8}, {6, -8}}, color = {0, 127, 255}));
      connect(heater.port_a, HTRecuperator1.from_comp_port_b) annotation(
        Line(points = {{32, 20}, {32, 20}, {32, -8}, {32, -8}}, color = {0, 127, 255}));
      connect(turbine.port_b, HTRecuperator1.from_turb_port_a) annotation(
        Line(points = {{80, 4}, {80, 4}, {80, -18}, {32, -18}, {32, -18}}, color = {0, 127, 255}));
      connect(heater.port_b, turbine.port_a) annotation(
        Line(points = {{32, 36}, {70, 36}, {70, 14}, {70, 14}}, color = {0, 127, 255}));
    
// Overall equations
      P_nom = (-turbine.W_turb) - firstCompressor.W_comp - reCompressor.W_comp - secondCompressor.W_comp;
      efficiencyCycle * heater.Q_heater = P_nom;
      E_bal_check = turbine.W_turb + firstCompressor.W_comp + reCompressor.W_comp + secondCompressor.W_comp + heater.Q_heater + cooler.Q_cooler + secondCooler.Q_cooler;
      W_out = P_nom / turbine.port_a.m_flow;
      UA_LTR = sum(LTRecuperator.UA_dis);
      UA_HTR = sum(HTRecuperator1.UA_dis);
// Exergy analysis
//  ex_d_tot = mainCompressor.ex_d + LTRecuperator.ex_d + HTRecuperator1.ex_d + reCompressor.ex_d + heater.ex_d + turbine.ex_d + cooler.ex_d;
//  mainCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_mainCompressor;
//  LTRecuperator.ex_d * 100 / ex_d_tot = ex_d_percent_LTRecuperator;
//  HTRecuperator1.ex_d * 100 / ex_d_tot = ex_d_percent_HTRecuperator;
//  reCompressor.ex_d * 100 / ex_d_tot = ex_d_percent_reCompressor;
//  heater.ex_d * 100 / ex_d_tot = ex_d_percent_heater;
//  turbine.ex_d * 100 / ex_d_tot = ex_d_percent_turbine;
//  cooler.ex_d * 100 / ex_d_tot = ex_d_percent_cooler;
//  ex_in = heater.Q_heater * (1 - T_amb / T_high);
//  eta_ex = 1 - ex_d_tot / ex_in;
//  eta_carnot = 1 - T_amb / T_high;
      annotation(
        Icon);
    end ModifiedRecompCycle;