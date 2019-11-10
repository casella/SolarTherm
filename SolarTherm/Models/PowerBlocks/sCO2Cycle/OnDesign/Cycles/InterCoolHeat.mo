within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cycles;

  model InterCoolHeat
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      //Parameters
      parameter Modelica.SIunits.AbsolutePressure p_high = 250 * 10 ^ 5 "high pressure of the cycle";
      parameter Modelica.SIunits.ThermodynamicTemperature T_high = 775 + 273.15 "inlet temperature of the turbine";
      parameter Modelica.SIunits.ThermodynamicTemperature T_amb = 40 + 273.15 "ambiant  temperature";
      parameter Real PR = 2.313 "Pressure ratio";
      parameter Modelica.SIunits.Power P_nom = 10 ^ 8 "Nominal electrical power";
      parameter Modelica.SIunits.Efficiency eta_comp = 0.87 "Isentropic efficiency of the compressors";
      parameter Modelica.SIunits.Efficiency eta_turb = 0.9 "Isentropic efficiency of the turbine";
      parameter Real gamma = 0.2 "Part of the mass flow going to the recompression directly";
      parameter Integer N_q = 15 "Number of discretization of the heat recuperators";
      // Variables to investigate the cycle and its simulation.
      Modelica.SIunits.Efficiency efficiencyCycle "Efficiency of the cycle";
      Real E_bal_check;
      SolarTherm.Types.SpecificWork W_out "Specific Work of the cycle";
      // Instanciation of the components
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft comprCool(PR = 1.6, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-28, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler interCooler(T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-50, -22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Cooler mainCooler(T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-74, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowSplitter splitter(gamma = gamma) annotation(
        Placement(visible = true, transformation(origin = {-2, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft comprRecomp(PR = PR, eta_comp = eta_comp, p_out = p_high, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {2, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve LTRecuperator(N_q = N_q, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {32, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.FlowMixer mixer annotation(
        Placement(visible = true, transformation(origin = {60, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.HeatRecuperatorDTAve HTRecuperator(N_q = N_q, flowGuess = P_nom / 10 ^ 5, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {88, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Heater mainHeater(T_high = T_high, T_amb = T_amb, is_first_heater = true) annotation(
        Placement(visible = true, transformation(origin = {34, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine mainTurbine(PR = PR / 1.6, eta_turb = eta_turb, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {38, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Turbine reHeatTurbine(PR = 1.6, eta_turb = eta_turb, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {74, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.Heater reHeater(T_high = T_high, T_amb = T_amb, is_first_heater=false) annotation(
        Placement(visible = true, transformation(origin = {58, -44}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign.CompressorOnShaft mainCompr(PR = PR / 1.6, eta_comp = eta_comp, p_out = p_high / PR * 1.6, T_amb = T_amb) annotation(
        Placement(visible = true, transformation(origin = {-74, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      efficiencyCycle * (mainHeater.Q_heater + reHeater.Q_heater) = P_nom;
      E_bal_check = mainTurbine.W_turb + reHeatTurbine.W_turb + mainCompr.W_comp + comprRecomp.W_comp + comprCool.W_comp + mainHeater.Q_heater + mainCooler.Q_cooler + reHeater.Q_heater + interCooler.Q_cooler;
      W_out = P_nom / mainTurbine.port_a.m_flow;
      P_nom = (-mainTurbine.W_turb) - reHeatTurbine.W_turb - mainCompr.W_comp - comprRecomp.W_comp - comprCool.W_comp;
      connect(mainCompr.port_b, interCooler.port_a) annotation(
        Line(points = {{-70, -44}, {-60, -44}, {-60, -22}, {-58, -22}, {-58, -22}}, color = {0, 127, 255}));
      connect(mainCooler.port_b, mainCompr.port_a) annotation(
        Line(points = {{-74, -10}, {-80, -10}, {-80, -32}, {-82, -32}, {-82, -32}}, color = {0, 127, 255}));
      connect(reHeatTurbine.port_b, HTRecuperator.from_turb_port_a) annotation(
        Line(points = {{80, -74}, {92, -74}, {92, -2}, {94, -2}, {94, -2}}, color = {0, 127, 255}));
      connect(reHeater.port_b, reHeatTurbine.port_a) annotation(
        Line(points = {{66, -44}, {72, -44}, {72, -64}, {70, -64}, {70, -64}}, color = {0, 127, 255}));
      connect(mainTurbine.port_b, reHeater.port_a) annotation(
        Line(points = {{44, -74}, {50, -74}, {50, -44}}, color = {0, 127, 255}));
      connect(mainHeater.port_b, mainTurbine.port_a) annotation(
        Line(points = {{34, -46}, {34, -46}, {34, -64}, {34, -64}}, color = {0, 127, 255}));
      connect(HTRecuperator.from_comp_port_b, mainHeater.port_a) annotation(
        Line(points = {{94, 8}, {34, 8}, {34, -30}}, color = {0, 127, 255}));
      connect(mixer.port_b, HTRecuperator.from_comp_port_a) annotation(
        Line(points = {{60, 4}, {80, 4}, {80, 8}, {80, 8}}, color = {0, 127, 255}));
      connect(LTRecuperator.from_comp_port_a, comprCool.port_b) annotation(
        Line(points = {{24, 26}, {-22, 26}, {-22, -50}, {-24, -50}, {-24, -50}}, color = {0, 127, 255}));
      connect(LTRecuperator.from_turb_port_b, splitter.port_a) annotation(
        Line(points = {{24, 16}, {6, 16}, {6, 2}, {6, 2}}, color = {0, 127, 255}));
      connect(HTRecuperator.from_turb_port_b, LTRecuperator.from_turb_port_a) annotation(
        Line(points = {{80, -2}, {72, -2}, {72, 16}, {38, 16}, {38, 16}}, color = {0, 127, 255}));
      connect(LTRecuperator.from_comp_port_b, mixer.second_port_a) annotation(
        Line(points = {{38, 26}, {52, 26}, {52, -4}, {52, -4}, {52, -4}}, color = {0, 127, 255}));
      connect(mixer.first_port_a, comprRecomp.port_b) annotation(
        Line(points = {{60, -12}, {60, -12}, {60, -24}, {4, -24}, {4, -52}, {6, -52}}, color = {0, 127, 255}));
      connect(splitter.one_gamma_port_b, mainCooler.port_a) annotation(
        Line(points = {{-10, 2}, {-60, 2}, {-60, 6}, {-74, 6}, {-74, 6}}, color = {0, 127, 255}));
      connect(splitter.gamma_port_b, comprRecomp.port_a) annotation(
        Line(points = {{-2, 10}, {-6, 10}, {-6, -40}, {-6, -40}}, color = {0, 127, 255}));
      connect(interCooler.port_b, comprCool.port_a) annotation(
        Line(points = {{-42, -22}, {-34, -22}, {-34, -38}, {-36, -38}, {-36, -38}}, color = {0, 127, 255}));
    
      annotation(
        Documentation(info = "<html>
    		<p> This is the model of an intercooling cycle with recompression and reheating. It is here to show how convenient it can be to just connect the elements. </p>
    		<p> Exergy and financial analysis are not implemented, they can be taken from another model. </p>
    		</html>"));
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end InterCoolHeat;