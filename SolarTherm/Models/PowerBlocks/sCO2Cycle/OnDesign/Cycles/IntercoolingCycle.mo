within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model IntercoolingCycle "Model of an intercooling cycle non optimized"
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      //Parameters
      parameter Real deltaT_cooler = 15 "Approach difference of temperature at the outlet of the cooler";
      parameter Modelica.SIunits.AbsolutePressure p_high = 200 * 10 ^ 5 "high pressure of the cycle";
      parameter Modelica.SIunits.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the turbine";
      parameter Modelica.SIunits.ThermodynamicTemperature T_amb = 40 + 273.15 "ambiant temperature";
      parameter Real PR = 2.371 "Pressure ratio";
      parameter Modelica.SIunits.Power P_nom = 10 ^ 8 "Nominal electrical power";
      parameter Modelica.SIunits.Efficiency eta_comp = 0.87 "Isentropic efficiency of the compressors";
      parameter Modelica.SIunits.Efficiency eta_turb = 0.9 "Isentropic efficiency of the turbine";
    
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
      // Instanciation of the components
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Heater heater(T_high = T_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {32, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine turbine(PR = PR, eta_turb = eta_turb, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {74, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft compressor1(PR = 1.5, eta_comp = eta_comp, p_out = p_high / PR * 1.5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-66, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft compressor2(PR = PR / 1.5, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-14, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler1(T_amb = T_amb, deltaT_cooler = deltaT_cooler) annotation(
        Placement(visible = true, transformation(origin = {-62, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler cooler2(T_amb = T_amb, deltaT_cooler = deltaT_cooler) annotation(
        Placement(visible = true, transformation(origin = {-40, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTRecuperator(N_q = N_q, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-2, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Variables to investigate the cycle and its simulation.
      Modelica.SIunits.Efficiency efficiencyCycle;
      Real E_bal_check;
      //Exergy analysis
      Real ex_d_percent_compressor1 "First compressor exergy destruction";
      Real ex_d_percent_compressor2 "Second compressor exergy destruction";
      Real ex_d_percent_HTRecuperator "Recuperator exergy destruction";
      Real ex_d_percent_heater "heater exergy destruction";
      Real ex_d_percent_turbine "turbine exergy destruction";
      Real ex_d_percent_cooler1 "First cooler exergy destruction";
      Real ex_d_percent_cooler2 "Second cooler exergy destruction";
      SolarTherm.Types.SpecificExergy ex_d_tot "Total exergy destruction";
      SolarTherm.Types.SpecificExergy ex_in;
      Modelica.SIunits.Efficiency eta_ex "Exergetic efficiency";
    equation
//Connectors
      connect(compressor2.port_b, HTRecuperator.from_comp_port_a) annotation(
        Line(points = {{-10, 14}, {-10, 14}, {-10, -4}, {-10, -4}}, color = {0, 127, 255}));
      connect(cooler2.port_b, compressor2.port_a) annotation(
        Line(points = {{-40, 30}, {-22, 30}, {-22, 26}, {-22, 26}, {-22, 26}}, color = {0, 127, 255}));
      connect(compressor1.port_b, cooler2.port_a) annotation(
        Line(points = {{-62, 14}, {-40, 14}, {-40, 14}, {-40, 14}}, color = {0, 127, 255}));
      connect(cooler1.port_b, compressor1.port_a) annotation(
        Line(points = {{-62, 0}, {-82, 0}, {-82, 26}, {-74, 26}, {-74, 26}, {-74, 26}}, color = {0, 127, 255}));
      connect(HTRecuperator.from_turb_port_b, cooler1.port_a) annotation(
        Line(points = {{-10, -14}, {-10, -14}, {-10, -16}, {-62, -16}, {-62, -16}}, color = {0, 127, 255}));
      connect(HTRecuperator.from_turb_port_a, turbine.port_b) annotation(
        Line(points = {{4, -14}, {80, -14}, {80, 4}}, color = {0, 127, 255}));
      connect(HTRecuperator.from_comp_port_b, heater.port_a) annotation(
        Line(points = {{4, -4}, {4, 7}, {32, 7}, {32, 20}}, color = {0, 127, 255}));
      connect(heater.port_b, turbine.port_a) annotation(
        Line(points = {{32, 36}, {70, 36}, {70, 14}, {70, 14}}, color = {0, 127, 255}));
// cycle efficiency; mass flow is imposed with the Power equation.
      P_nom = (-turbine.W_turb) - compressor1.W_comp - compressor2.W_comp;
      efficiencyCycle * heater.Q_heater = P_nom;
      E_bal_check = turbine.W_turb + compressor1.W_comp + compressor2.W_comp + heater.Q_heater + cooler1.Q_cooler + cooler2.Q_cooler;
// Exergy analysis
      ex_d_tot = compressor1.ex_d + compressor2.ex_d + HTRecuperator.ex_d + heater.ex_d + turbine.ex_d + cooler1.ex_d + cooler2.ex_d;
      compressor1.ex_d * 100 / ex_d_tot = ex_d_percent_compressor1;
      compressor2.ex_d * 100 / ex_d_tot = ex_d_percent_compressor2;
      HTRecuperator.ex_d * 100 / ex_d_tot = ex_d_percent_HTRecuperator;
      heater.ex_d * 100 / ex_d_tot = ex_d_percent_heater;
      turbine.ex_d * 100 / ex_d_tot = ex_d_percent_turbine;
      cooler1.ex_d * 100 / ex_d_tot = ex_d_percent_cooler1;
      cooler2.ex_d * 100 / ex_d_tot = ex_d_percent_cooler2;
      ex_in = heater.Q_heater * (1 - T_amb / T_high);
      eta_ex = ((-turbine.W_turb) - compressor1.W_comp - compressor2.W_comp) / ex_in;
    end IntercoolingCycle;