within SolarTherm.Models.PowerBlocks.sCO2Cycle;

package QuickCO2PB

  model recompPB
    extends SolarTherm.Media.CO2.PropCO2;
    import SI = Modelica.SIunits;
    import FI = SolarTherm.Models.Analysis.Finances;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    //input SolarTherm.Interfaces.Connectors.WeatherBus wbus;
    extends Icons.PowerBlock;
    //  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
    //    Placement(visible = true,transformation(extent = {{-54, 22}, {-34, 42}}, rotation = 0), iconTransformation(extent = {{-48, 28}, {-40, 36}}, rotation = 0)));
    //  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
    //    Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
    Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
      Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{48, -12}, {58, -2}}, rotation = 0)));
    // PB parameters
    parameter Boolean external_parasities = false "= true enable parasities as an input";
    parameter Real nu_min = 0.25 "Minimum turbine operation";
    Modelica.Blocks.Interfaces.RealInput parasities if external_parasities annotation(
      Placement(transformation(extent = {{-12, -12}, {12, 12}}, rotation = -90, origin = {1.77636e-015, 80}), iconTransformation(extent = {{-6, -6}, {6, 6}}, rotation = -90, origin = {20, 60})));
    input SI.ThermodynamicTemperature T_amb;
    //Cycle parameters
    parameter SI.AbsolutePressure p_high = 200 * 10 ^ 5 "high pressure of the cycle";
    parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the turbine";
    parameter SI.ThermodynamicTemperature T_amb_des = 30 + 273.15 "ambiant temperature";
    parameter Real PR = 2.5 "Pressure ratio";
    parameter SI.Power P_gro = 100 * 10 ^ 6 "first guess of power outlet";
    parameter SI.Power P_nom(fixed = false) "Electrical power at design point";
    parameter SI.MassFlowRate m_HTF_des = 1000 "Mass flow rate at design point";
    parameter Real gamma = 0.28 "Part of the mass flow going to the recompression directly";
    parameter SI.AngularVelocity[4] choiceN = {75000, 30000, 10000, 3600} * 0.10471975512;
    parameter SI.AngularVelocity N_shaft = choiceN[integer(Modelica.Math.log(P_gro / 10 ^ 6) / Modelica.Math.log(10)) + 2];
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
    parameter Real ratio_m_des = 1 - gamma;
    //Cooler parameters
    parameter SI.ThermodynamicTemperature T_low = 45 + 273.15 "Outlet temperature of the cooler";
    //Exchanger parameters
    parameter SI.ThermodynamicTemperature T_HTF_in_des = 800 + 273.15;
    parameter Integer N_exch = 5;
    //Financial analysis
    parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
    parameter FI.Money C_LTR(fixed = false) "cost of the low temperature heat recuperator";
    parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
    parameter FI.Money C_mainCompressor(fixed = false) "cost of the main compressor";
    parameter FI.Money C_reCompressor(fixed = false) "cost of the re compressor";
    parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
    parameter FI.Money C_generator(fixed = false) "cost of the generator";
    parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
    parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
    parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
    //Results
    SI.Efficiency eta_cycle;
    SI.Energy E_net(final start = 0, fixed = true, displayUnit = "MW.h");
    Boolean m_sup "Disconnect the production of electricity when the outlet pressure of the turbine is close to the critical pressure";
    //Components instanciation
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve HTR(N_q = N_HTR, P_nom_des = P_gro, ratio_m_des = 1) annotation(
      Placement(visible = true, transformation(origin = {12, -22}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.CompressorOnShaft mainCompressor(eta_design = eta_comp_main, N_design = N_shaft, P_nom_des = P_gro, p_high_des = p_high) annotation(
      Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Cooler cooler(T_low = T_low, P_nom_des = P_gro, T_amb_des = T_amb_des) annotation(
      Placement(visible = true, transformation(origin = {-78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Turbine turbine(PR = PR, N_shaft = N_shaft, eta_design = eta_turb) annotation(
      Placement(visible = true, transformation(origin = {66, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Exchanger exchanger(redeclare package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph, P_nom_des = P_gro, T_out_CO2_des = T_high, N_exch = N_exch, ratio_m_des = 1) annotation(
      Placement(visible = true, transformation(extent = {{34, -8}, {54, 12}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.CompressorOnShaft reCompressor(N_design = N_shaft, P_nom_des = P_gro, p_high_des = p_high) annotation(
      Placement(visible = true, transformation(origin = {-54, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve LTR(N_q = N_LTR, P_nom_des = P_gro, ratio_m_des = 1 - gamma) annotation(
      Placement(visible = true, transformation(origin = {-42, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.FlowMixer mixer annotation(
      Placement(visible = true, transformation(origin = {-20, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.FlowSplitter splitter(gamma = gamma) annotation(
      Placement(visible = true, transformation(origin = {-62, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter MedRec.ThermodynamicState state_HTF_in_des = MedRec.setState_pTX(1.0325 * 10 ^ 5, T_HTF_in_des);
    //     SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = 800+273.15, p_out = 10 ^ 5, m_flow = exchanger.m_HTF_des, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
    //        Placement(visible = true, transformation(origin = {-52, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    //      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sink annotation(
    //        Placement(visible = true, transformation(origin = {-50, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // Modelica.Blocks.Interfaces.RealInput parasities_internal;
  protected
  initial equation
    exchanger.h_in_HTF_des = MedRec.specificEnthalpy(state_HTF_in_des);
    exchanger.p_in_HTF_des = state_HTF_in_des.p;
    exchanger.m_HTF_des = m_HTF_des;
    P_nom = (-turbine.W_turb_des) - mainCompressor.W_comp_des - reCompressor.W_comp_des - cooler.P_cool_des;
// enthalpy equalities
//main loop
    exchanger.h_in_CO2_des = HTR.h_out_comp_des;
    turbine.h_in_des = exchanger.h_out_CO2_des;
    HTR.h_in_turb_des = turbine.h_out_des;
    LTR.h_in_turb_des = HTR.h_out_turb_des;
    cooler.h_in_des = LTR.h_out_turb_des;
    mainCompressor.h_in_des = cooler.h_out_des;
    LTR.h_in_comp_des = mainCompressor.h_out_des;
// recompression loop
    reCompressor.h_in_des = LTR.h_out_turb_des;
    HTR.h_in_comp_des = ratio_m_des * LTR.h_out_comp_des + (1 - ratio_m_des) * reCompressor.h_out_des;
//pressure equalities
//main loop
    exchanger.p_in_CO2_des = HTR.p_out_comp_des;
    turbine.p_in_des = exchanger.p_out_CO2_des;
    HTR.p_in_turb_des = turbine.p_out_des;
    LTR.p_in_turb_des = HTR.p_out_turb_des;
    cooler.p_in_des = LTR.p_out_turb_des;
    mainCompressor.p_in_des = cooler.p_out_des;
    LTR.p_in_comp_des = mainCompressor.p_out_des;
//recompression loop
    reCompressor.p_in_des = LTR.p_out_turb_des;
    HTR.p_in_comp_des = ratio_m_des * LTR.p_out_comp_des + (1 - ratio_m_des) * reCompressor.p_out_des;
//mass flow equalities
//main loop
//exchanger.m_CO2_des = HTR.m_comp_des;
    turbine.m_des = exchanger.m_CO2_des;
    HTR.m_turb_des = turbine.m_des;
    LTR.m_turb_des = HTR.m_turb_des;
    cooler.m_des = LTR.m_turb_des * ratio_m_des;
    mainCompressor.m_des = cooler.m_des;
    LTR.m_comp_des = mainCompressor.m_des;
//recompression loop
    HTR.m_comp_des = reCompressor.m_des + LTR.m_comp_des;
    reCompressor.m_des = gamma * LTR.m_turb_des;
// Financial Analysis
    C_HTR = if HTR.T_turb_des[N_HTR] >= 550 + 273.15 then 49.45 * HTR.UA_HTR ^ 0.7544 * (1 + 0.02141 * (HTR.T_turb_des[N_HTR] - 550 - 273.15)) else 49.45 * HTR.UA_HTR ^ 0.7544;
    C_LTR = 49.45 * LTR.UA_HTR ^ 0.7544;
    C_turbine = if exchanger.T_CO2_des[2] >= 550 + 273.15 then 406200 * (-turbine.W_turb_des / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (exchanger.T_CO2_des[2] - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb_des / 10 ^ 6) ^ 0.8;
    C_mainCompressor = 1230000 * (mainCompressor.W_comp_des / 10 ^ 6) ^ 0.3392;
    C_reCompressor = 1230000 * (reCompressor.W_comp_des / 10 ^ 6) ^ 0.3392;
    C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
    C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
    C_exchanger = pri_exchanger * exchanger.Q_HX_des * m_HTF_des / 1000;
    C_PB = (C_HTR + C_LTR + C_turbine + C_mainCompressor + C_reCompressor + C_generator + C_cooler + C_exchanger) * 1.05;
// 1.05 corresponds to inflation from 2017, as correlations are in 2017' dollars.
  equation
    connect(m_sup, exchanger.m_sup) annotation(
      Line);
    connect(exchanger.CO2_port_b, turbine.port_a) annotation(
      Line(points = {{50, -4}, {64.2, -4}, {64.2, -23.8}, {62.2, -23.8}, {62.2, -25.8}}, color = {0, 127, 255}));
    connect(exchanger.CO2_port_a, HTR.from_comp_port_b) annotation(
      Line(points = {{37, -4}, {22.8, -4}, {22.8, -13.6}}, color = {0, 127, 255}));
//  connect(fluid_b, exchanger.HTF_port_b) annotation(
//    Line(points = {{-64, -50}, {12, -50}, {12, 10}, {38, 10}, {38, 10}}, color = {0, 127, 255}));
//  connect(fluid_a, exchanger.HTF_port_a) annotation(
//    Line(points = {{-44, 32}, {54, 32}, {54, 10}, {52, 10}, {52, 10}}, color = {0, 127, 255}));
//connect(src.port_b, fluid_a) annotation(
//  Line(points = {{-60, 64}, {-72, 64}, {-72, 30}, {-42, 30}, {-42, 32}, {-44, 32}}, color = {0, 127, 255}));
//connect(sink.port_a, fluid_b) annotation(
//  Line(points = {{-42, -78}, {-64, -78}, {-64, -54}, {-64, -54}, {-64, -50}}, color = {0, 127, 255}));
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
    connect(turbine.port_b, HTR.from_turb_port_a) annotation(
      Line(points = {{72, -36.6}, {48, -36.6}, {48, -34.6}, {22, -34.6}, {22, -34.1}, {22, -34.1}, {22, -31.6}}, color = {0, 127, 255}));
    connect(cooler.port_b, mainCompressor.port_a) annotation(
      Line(points = {{-78, -46}, {-88, -46}, {-88, -6}, {-82, -6}, {-82, -6}}, color = {0, 127, 255}));
    connect(cooler.T_amb, T_amb);
    connect(m_sup, cooler.m_sup);
//  exchanger.HTF_port_b.m_flow=fluid_b.m_flow;
//  -fluid_a.m_flow+exchanger.HTF_port_a.m_flow=0;
//  exchanger.HTF_port_b.p=fluid_b.p;
//  exchanger.HTF_port_a.p=fluid_a.p;
//  exchanger.HTF_port_a.h_outflow=fluid_a.h_outflow;
//  fluid_b.h_outflow=exchanger.HTF_port_b.h_outflow;
//m_sup = true;
    m_sup = exchanger.HTF_port_a.m_flow >= exchanger.m_HTF_des * nu_min;
    exchanger.CO2_port_a.m_flow = exchanger.m_CO2_des;
    eta_cycle = W_net / exchanger.Q_HX;
    der(E_net) = W_net;
    W_net = if m_sup then (-turbine.W_turb) - mainCompressor.W_comp - reCompressor.W_comp - cooler.P_cooling else 0;
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end recompPB;

  model test
    extends SolarTherm.Media.CO2.PropCO2;
    import CV = Modelica.SIunits.Conversions;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    ComplicatedPB.CompressorOnShaft compr annotation(
      Placement(visible = true, transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation

    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end test;

  model interpolSimplePB "This model calculates a few off-design points as initial equation and after interpolate over them"
    import SI = Modelica.SIunits;
    import CV = Modelica.SIunits.Conversions;
    import FI = SolarTherm.Models.Analysis.Finances;
    extends SolarTherm.Media.CO2.PropCO2;
    extends Icons.PowerBlock;
    input SI.ThermodynamicTemperature T_amb;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    //Objective: having a few points in order to interpolate the HTF outlet temperature and the TIT. Should be feasible with for loops in the initial equation.
    parameter SI.Power P_nom = 100 * 10 ^ 6;
    parameter SI.AbsolutePressure p_high = CV.from_bar(250);
    parameter SI.AbsolutePressure p_low = CV.from_bar(80);
    parameter Real PR = p_high / p_low;
    parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
    parameter SI.ThermodynamicTemperature T_HTF_in_des = CV.from_degC(800) "inlet temperature of the HTF";
    parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
    parameter SI.ThermodynamicTemperature T_amb_des = CV.from_degC(20) "ambient temperature at design";
    parameter SI.TemperatureDifference pinchHTR = 5;
    parameter SI.TemperatureDifference pinchExch = T_HTF_in_des - T_high;
    parameter SI.Efficiency eta_comp = 0.89;
    parameter SI.Efficiency eta_turb = 0.93;
    parameter Integer N_exch = 5;
    parameter Integer N_HTR = 15;
    parameter Integer N_cooler = 15;
    //Variables for real time calculation
    //Exchanger variables
    Boolean m_sup "indicates if mass flow from tank is superior enough to calculate";
    SI.MassFlowRate m_HTF_bis "mass flow rate used for switching off the PB";
    SI.ThermodynamicTemperature T_HTF_in "Inlet temperature of HTF";
    SI.ThermodynamicTemperature T_HTF_out "Outlet temperature of HTF";
    SI.ThermodynamicTemperature TIT_CO2 "Turbine inlet Temperature of CO2";
    SI.ThermodynamicTemperature T_CO2_in_HX "CO2 temperature at the inlet of the exchanger";
    MedPB.ThermodynamicState[N_exch] state_CO2_HX(each h.start = 10 ^ 6);
    MedRec.ThermodynamicState[N_exch] state_HTF_HX(each h.start = 6 * 10 ^ 5);
    SI.TemperatureDifference[N_exch] deltaT_HX(each start = 20);
    //PowerBlock variables
    Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
      Placement(transformation(extent = {{-54, 22}, {-34, 42}}), iconTransformation(extent = {{-48, 30}, {-40, 38}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
      Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
    Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
      Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
    SI.Power P_cooling "Cooling power necessary to cool down the fluid";
    parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
    parameter SI.ThermodynamicTemperature T_HTF_out_des(fixed = false);
    parameter MedPB.ThermodynamicState state_comp_out(p.fixed = false, h.fixed = false);
    //Off-design parameters for calculation
    parameter Integer N_points = 10;
    parameter SI.ThermodynamicTemperature[N_points] T_high_frac(each fixed = false) "points for which efficiency will be calculated";
    //parameters to calculate on and off-design
    parameter SI.MassFlowRate m_des(fixed = false, start = P_nom / 10 ^ 5) "mass flow rate at the outlet of the turbine";
    parameter SI.Power[N_points] W_out(each fixed = false, each start = P_nom) "power at the outlet for off-design";
    parameter SI.Efficiency eta_cycle_on(fixed = false, start = 0.5) "efficiency of the cycle";
    parameter SI.Efficiency[N_points] eta_cycle_off(each fixed = false, each start = 0.5) "efficiency of the cycle";
    // bougé pour voir ce qu'il se passe
    parameter SI.ThermodynamicTemperature[N_points] T_CO2_in_HX_table(each fixed = false);
    parameter SI.Power W_comp(fixed = false);
    parameter SI.HeatFlowRate Q_dis_exch(fixed = false);
    parameter MedPB.ThermodynamicState state_cooler_out(p.fixed = false, h.fixed = false);
  protected
    //Inlet of the turbine
    parameter SI.SpecificEntropy[N_points + 1] s_in_turb(each fixed = false, each start = 2900) "inlet entropy to the turbine";
    parameter MedPB.ThermodynamicState[N_points + 1] state_in_turb(each p.fixed = false, each h.fixed = false);
    parameter SI.Power[N_points + 1] W_turb(each fixed = false);
    parameter MedPB.ThermodynamicState[N_points + 1] state_out_turb(each p.fixed = false, each h.fixed = false);
    //Main Compressor parameters
    parameter SI.SpecificEntropy s_in_mainComp(fixed = false, start = 1500) "inlet entropy to the main compressor";
    //Exchanger parameters. Calculated only at on-design
    parameter SolarTherm.Types.Conductance[N_exch - 1] UA_exch_dis(each fixed = false, each start = P_nom * 0.4 / N_exch);
    parameter MedPB.ThermodynamicState[N_exch] state_CO2_exch(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
    parameter MedRec.ThermodynamicState[N_exch] state_HTF_exch(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
    parameter SI.TemperatureDifference[N_exch] deltaT_exch(each start = 20, each fixed = false);
    // HTR parameters
    parameter SolarTherm.Types.Conductance[N_HTR - 1] UA_HTR_dis(each fixed = false, each start = P_nom * 0.4 / N_HTR);
    parameter SI.HeatFlowRate Q_dis_HTR(fixed = false, start = 2 * 10 ^ 4);
    parameter MedPB.ThermodynamicState[N_points + 1, N_HTR] state_HTR_turb(each p.fixed = false, each h.fixed = false, each h.start = 7 * 10 ^ 5);
    parameter MedPB.ThermodynamicState[N_points + 1, N_HTR] state_HTR_comp(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
    parameter SI.TemperatureDifference[N_points + 1, N_HTR] deltaT_HTR(each start = 20, each fixed = false);
    // Cooler parameters
    parameter MedPB.ThermodynamicState[N_cooler] state_cooler(each p.fixed = false, each h.fixed = false, each h.start = 4.5 * 10 ^ 5);
    parameter SI.Power P_cool_des(fixed = false) "on-design power necessary to run the fans";
    parameter SI.HeatFlowRate Q_dis_cooler(fixed = false);
    parameter SolarTherm.Types.Conductance[N_cooler - 1] UA_cooler_dis(each fixed = false);
    //Financial analysis
    parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
    parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
    parameter FI.Money C_compressor(fixed = false) "cost of the main compressor";
    parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
    parameter FI.Money C_generator(fixed = false) "cost of the generator";
    parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
    parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
  initial equation
//Part I. On-design=first point
//Inlet of the turbine
    state_in_turb[1] = MedPB.setState_pTX(p_high, T_high);
    s_in_turb[1] = MedPB.specificEntropy(state_in_turb[1]);
    state_out_turb[1] = MedPB.setState_phX(p_low, state_in_turb[1].h + eta_turb * (MedPB.specificEnthalpy(MedPB.setState_psX(p_low, s_in_turb[1])) - state_in_turb[1].h));
    W_turb[1] = m_des * (state_out_turb[1].h - state_in_turb[1].h);
//High temperature recuperator
    for k in 1:N_HTR - 1 loop
      Q_dis_HTR = state_HTR_turb[1, k + 1].h - state_HTR_turb[1, k].h;
      Q_dis_HTR = state_HTR_comp[1, k + 1].h - state_HTR_comp[1, k].h;
      deltaT_HTR[1, k] = MedPB.temperature(state_HTR_turb[1, k]) - MedPB.temperature(state_HTR_comp[1, k]);
      m_des * Q_dis_HTR = UA_HTR_dis[k] * (deltaT_HTR[1, k] + deltaT_HTR[1, k + 1]) / 2;
      state_HTR_comp[1, k + 1].p = p_high;
//need to use setState when doing isentropic enthalpy
      state_HTR_turb[1, k].p = p_low;
    end for;
    deltaT_HTR[1, N_HTR] = MedPB.temperature(state_HTR_turb[1, N_HTR]) - MedPB.temperature(state_HTR_comp[1, N_HTR]);
    min(deltaT_HTR[1]) = pinchHTR;
    state_HTR_turb[1, N_HTR].p = p_low;
    state_HTR_turb[1, N_HTR].h = state_out_turb[1].h;
    state_cooler[1].h = state_HTR_turb[1, 1].h;
//Cooler. UA (and therefore this part) are only for economical correlation
    for k in 1:N_cooler - 1 loop
      Q_dis_cooler = state_cooler[k + 1].h - state_cooler[k].h;
      m_des * Q_dis_cooler = -UA_cooler_dis[k] * (MedPB.temperature(state_cooler[k + 1]) - T_amb_des + MedPB.temperature(state_cooler[k]) - T_amb_des) / 2;
      state_cooler[k].p = p_low;
    end for;
    state_cooler[N_cooler] = state_cooler_out;
    state_cooler_out = MedPB.setState_pTX(p_low, T_low);
    P_cool_des * (T_low - T_amb_des) / (-m_des * Q_dis_cooler * (N_cooler - 1)) = 1.49 * 10 ^ 6 * (35.7 - 30) / (136.6 * 10 ^ 6);
//Main compressor. Inlet supposed always the same for off-design
    s_in_mainComp = MedPB.specificEntropy(state_cooler_out);
    state_HTR_comp[1, 1] = MedPB.setState_phX(p_high, state_cooler_out.h + (MedPB.specificEnthalpy(MedPB.setState_psX(p_high, s_in_mainComp)) - state_cooler_out.h) / eta_comp);
    W_comp = m_des * (state_HTR_comp[1, 1].h - state_cooler_out.h);
    state_comp_out = state_HTR_comp[1, 1];
//pour voir ce qu'il se passe
//Exchanger. Done only on-design. After, only the TIT is varied.
    state_CO2_exch[N_exch].h = state_in_turb[1].h;
    state_CO2_exch[1].h = state_HTR_comp[1, N_HTR].h;
    state_HTF_exch[N_exch] = MedRec.setState_pTX(10 ^ 5, T_HTF_in_des);
    for k in 1:N_exch - 1 loop
      Q_dis_exch = state_CO2_exch[k + 1].h - state_CO2_exch[k].h;
      Q_dis_exch = state_HTF_exch[k + 1].h - state_HTF_exch[k].h;
      m_des * Q_dis_exch = UA_exch_dis[k] * (deltaT_exch[k] + deltaT_exch[k + 1]) / 2;
      deltaT_exch[k] = MedRec.temperature(state_HTF_exch[k]) - MedPB.temperature(state_CO2_exch[k]);
      state_CO2_exch[k].p = p_high;
      state_HTF_exch[k].p = 10 ^ 5;
    end for;
    deltaT_exch[N_exch] = MedRec.temperature(state_HTF_exch[N_exch]) - MedPB.temperature(state_CO2_exch[N_exch]);
    T_HTF_out_des = MedRec.temperature(state_HTF_exch[1]);
    state_CO2_exch[N_exch].p = p_high;
    P_nom = (-W_turb[1]) - W_comp;
    eta_cycle_on = P_nom / (m_des * Q_dis_exch * (N_exch - 1));
//Part I.2. Financial analysis
    C_HTR = if MedPB.temperature(state_HTR_turb[1, N_exch]) >= 550 + 273.15 then 49.45 * sum(UA_HTR_dis) ^ 0.7544 * (1 + 0.02141 * (MedPB.temperature(state_HTR_turb[1, N_exch]) - 550 - 273.15)) else 49.45 * sum(UA_HTR_dis) ^ 0.7544;
    C_turbine = if T_high >= 550 + 273.15 then 406200 * (-W_turb[1] / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (T_high - 550 - 273.15) ^ 2) else 406200 * (-W_turb[1] / 10 ^ 6) ^ 0.8;
    C_compressor = 1230000 * (W_comp / 10 ^ 6) ^ 0.3392;
    C_cooler = 32.88 * sum(UA_cooler_dis) ^ 0.75;
    C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
    C_exchanger = pri_exchanger * m_des * Q_dis_exch * (N_exch - 1) / 1000;
    C_PB = (C_HTR + C_turbine + C_compressor + C_generator + C_cooler + C_exchanger) * 1.05;
//Part II. Off-design points. Mass flow rate is kept constant
    for i in 2:N_points + 1 loop
//    //Inlet of the turbine
      T_high_frac[i - 1] = (0.8 + (1.2 - 0.8) * (i - 1) / (N_points - 1)) * (T_high - 273.15) + 273.15;
      state_in_turb[i] = MedPB.setState_pTX(p_high, T_high_frac[i - 1]);
      s_in_turb[i] = MedPB.specificEntropy(state_in_turb[i]);
      state_out_turb[i] = MedPB.setState_phX(p_low, state_in_turb[i].h + eta_turb * (MedPB.specificEnthalpy(MedPB.setState_psX(p_low, s_in_turb[i])) - state_in_turb[i].h));
      state_HTR_turb[i, N_HTR] = state_out_turb[i];
      W_turb[i] = m_des * (state_out_turb[i].h - state_in_turb[i].h);
//High temperature recuperator
      for k in 1:N_HTR - 1 loop
        state_HTR_turb[i, k + 1].h - state_HTR_turb[i, k].h = state_HTR_comp[i, k + 1].h - state_HTR_comp[i, k].h;
        deltaT_HTR[i, k] = MedPB.temperature(state_HTR_turb[i, k]) - MedPB.temperature(state_HTR_comp[i, k]);
        m_des * (state_HTR_turb[i, k + 1].h - state_HTR_turb[i, k].h) = UA_HTR_dis[k] * (deltaT_HTR[i, k] + deltaT_HTR[i, k + 1]) / 2;
        state_HTR_turb[i, k].p = p_low;
        state_HTR_comp[i, k + 1].p = p_high;
      end for;
      deltaT_HTR[i, N_HTR] = MedPB.temperature(state_HTR_turb[i, N_HTR]) - MedPB.temperature(state_HTR_comp[i, N_HTR]);
      T_CO2_in_HX_table[i - 1] = MedPB.temperature(state_HTR_comp[i, N_HTR]);
//   //Inlet of the HTR. Main Compressor always has the same inlet so no off-design
      state_HTR_comp[i, 1] = MedPB.setState_phX(p_high, state_cooler[N_cooler].h + (MedPB.specificEnthalpy(MedPB.setState_psX(p_high, s_in_mainComp)) - state_cooler[N_cooler].h) / eta_comp);
      W_out[i - 1] = (-W_turb[i]) - W_comp;
      eta_cycle_off[i - 1] = W_out[i - 1] / (m_des * (state_in_turb[i].h - state_HTR_comp[i, N_HTR].h));
    end for;
  equation
//Exchanger code
    m_sup = fluid_a.m_flow >= 0.5 * m_des;
    m_HTF_bis = if m_sup then fluid_a.m_flow else m_des;
    for k in 1:N_exch - 1 loop
      m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = m_HTF_bis * (state_HTF_HX[k + 1].h - state_HTF_HX[k].h);
      m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = UA_exch_dis[k] * (deltaT_HX[k] + deltaT_HX[k + 1]) / 2;
      deltaT_HX[k] = MedRec.temperature(state_HTF_HX[k]) - MedPB.temperature(state_CO2_HX[k]);
      state_HTF_HX[k].p = state_HTF_HX[k + 1].p;
      state_CO2_HX[k].p = state_CO2_HX[k + 1].p;
    end for;
    deltaT_HX[N_exch] = MedRec.temperature(state_HTF_HX[N_exch]) - MedPB.temperature(state_CO2_HX[N_exch]);
    state_HTF_HX[N_exch] = if m_sup then MedRec.setState_phX(fluid_a.p, inStream(fluid_a.h_outflow)) else state_HTF_exch[N_exch];
//The TIT and T_CO2 at the inlet of the HX are interpolated from previous points
    T_HTF_in = MedRec.temperature(state_HTF_HX[N_exch]);
    T_HTF_out = MedRec.temperature(state_HTF_HX[1]);
    TIT_CO2 = MedPB.temperature(state_CO2_HX[N_exch]);
    T_CO2_in_HX = Modelica.Math.Vectors.interpolate(T_high_frac, T_CO2_in_HX_table, TIT_CO2);
    state_CO2_HX[1] = MedPB.setState_pTX(p_high, T_CO2_in_HX);
// Power Block calculations
    P_cooling = P_cool_des * ((T_low - T_amb_des) / (max(T_amb + 5, T_low) - T_amb)) ^ (3 / 0.805);
    W_net = if m_sup then max(0, Modelica.Math.Vectors.interpolate(T_high_frac, W_out, TIT_CO2) - P_cooling) else 0;
//Connectors obligations
    fluid_b.p = fluid_a.p;
    fluid_a.m_flow + fluid_b.m_flow = 0;
    fluid_b.h_outflow = state_HTF_HX[1].h;
    fluid_a.h_outflow = 0;
//shouldn't flow back
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end interpolSimplePB;

  model interpolWithCompo "This model calculates a few off-design points as initial equation and after interpolate over them"
    import SI = Modelica.SIunits;
    import CV = Modelica.SIunits.Conversions;
    import FI = SolarTherm.Models.Analysis.Finances;
    extends SolarTherm.Media.CO2.PropCO2;
    extends Icons.PowerBlock;
    input SI.ThermodynamicTemperature T_amb;
    replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
    replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
    //Objective: having a few points in order to interpolate the HTF outlet temperature and the TIT. Should be feasible with for loops in the initial equation.
    parameter SI.Power P_nom = 100 * 10 ^ 6;
    parameter SI.AbsolutePressure p_high = CV.from_bar(250);
    parameter SI.AbsolutePressure p_low = CV.from_bar(80);
    parameter Real PR = p_high / p_low;
    parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
    parameter SI.ThermodynamicTemperature T_HTF_in_des = CV.from_degC(800) "inlet temperature of the HTF";
    parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
    parameter SI.ThermodynamicTemperature T_amb_des = CV.from_degC(20) "ambient temperature at design";
    parameter SI.MassFlowRate m_des(fixed = false, start = P_nom / 10 ^ 5);
    parameter SI.TemperatureDifference pinchHTR = 5;
    parameter SI.TemperatureDifference pinchExch = T_HTF_in_des - T_high;
    parameter SI.Efficiency eta_comp = 0.89;
    parameter SI.Efficiency eta_turb = 0.93;
    parameter Integer N_exch = 5;
    parameter Integer N_HTR = 15;
    parameter Integer N_cooler = 15;
    //Variables for real time calculation
    //Exchanger variables
    Boolean m_sup "indicates if mass flow from tank is superior enough to calculate";
    SI.MassFlowRate m_HTF_bis "mass flow rate used for switching off the PB";
    SI.ThermodynamicTemperature T_HTF_in "Inlet temperature of HTF";
    SI.ThermodynamicTemperature T_HTF_out "Outlet temperature of HTF";
    SI.ThermodynamicTemperature TIT_CO2 "Turbine inlet Temperature of CO2";
    SI.ThermodynamicTemperature T_CO2_in_HX "CO2 temperature at the inlet of the exchanger";
    MedPB.ThermodynamicState[N_exch] state_CO2_HX(each h.start = 10 ^ 6);
    MedRec.ThermodynamicState[N_exch] state_HTF_HX(each h.start = 6 * 10 ^ 5);
    SI.TemperatureDifference[N_exch] deltaT_HX(each start = 20);
    //PowerBlock variables
    Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
      Placement(transformation(extent = {{-54, 22}, {-34, 42}}), iconTransformation(extent = {{-48, 30}, {-40, 38}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
      Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
    Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
      Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
    parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
    //Exchanger parameters. Calculated only at on-design
    parameter SolarTherm.Types.Conductance[N_exch - 1] UA_exch_dis(each fixed = false, each start = P_nom * 0.4 / N_exch);
    parameter SI.HeatFlowRate Q_dis_exch(fixed = false);
    parameter MedPB.ThermodynamicState[N_exch] state_CO2_exch(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
    parameter MedRec.ThermodynamicState[N_exch] state_HTF_exch(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
    parameter SI.TemperatureDifference[N_exch] deltaT_exch(each start = 20, each fixed = false);
    parameter SI.ThermodynamicTemperature T_HTF_out_des(fixed = false);
    //Off-design parameters for calculation
    parameter Integer N_points = 10;
    parameter SI.ThermodynamicTemperature[N_points] T_high_frac(each fixed = false) "TIT points for off-design";
    parameter SI.ThermodynamicTemperature[N_points] T_CO2_in_HX_table(each fixed = false) "Temperature at the inlet of the exchanger, for a TIT given. Used for interpolation";
    parameter SI.Power[N_points] W_out(each fixed = false) "Outlet power at off-design";
    parameter SI.Efficiency eta_cycle_on(fixed = false, start = 0.5) "efficiency of the cycle";
    parameter SI.Efficiency[N_points] eta_cycle_off(each fixed = false, each start = 0.5) "efficiency of the cycle";
    //Financial analysis
    parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
    parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
    parameter FI.Money C_compressor(fixed = false) "cost of the main compressor";
    parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
    parameter FI.Money C_generator(fixed = false) "cost of the generator";
    parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
    parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Turbine turbine(eta_design = eta_turb, PR = PR, N_points = N_points) annotation(
      Placement(visible = true, transformation(origin = {38, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Cooler cooler annotation(
      Placement(visible = true, transformation(origin = {-52, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve HTR(N_points = N_points, N_q = N_HTR, pinchRecuperator = pinchHTR, ratio_m_des = 1, P_nom = P_nom) annotation(
      Placement(visible = true, transformation(origin = {-2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    CompressorOnShaft compr(eta_design = eta_comp, PR = PR, N_points = N_points) annotation(
      Placement(visible = true, transformation(origin = {-38, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
//Part I. On-design=first point
//Inlet of the turbine
    turbine.state_in[1] = MedPB.setState_pTX(p_high, T_high);
    m_des = turbine.m_des;
    turbine.m_des = HTR.m_des;
//High temperature recuperator
    HTR.state_turb[1, N_HTR] = turbine.state_out[1];
    cooler.state_cooler[1] = HTR.state_turb[1, 1];
    cooler.m_des = HTR.m_des;
//Cooler. UA (and therefore this part) are only for economical correlation
    cooler.state_out = MedPB.setState_pTX(p_low, T_low);
//Main compressor. Inlet supposed always the same for off-design
    compr.state_in[1] = cooler.state_out;
    compr.m_des = cooler.m_des;
    HTR.state_comp[1, 1] = compr.state_out[1];
//Exchanger. Done only on-design. After, only the TIT is varied.
    state_CO2_exch[N_exch].h = turbine.state_in[1].h;
    state_CO2_exch[1].h = HTR.state_comp[1, N_HTR].h;
    state_HTF_exch[N_exch] = MedRec.setState_pTX(10 ^ 5, T_HTF_in_des);
    for k in 1:N_exch - 1 loop
      Q_dis_exch = state_CO2_exch[k + 1].h - state_CO2_exch[k].h;
      Q_dis_exch = state_HTF_exch[k + 1].h - state_HTF_exch[k].h;
      m_des * Q_dis_exch = UA_exch_dis[k] * (deltaT_exch[k] + deltaT_exch[k + 1]) / 2;
      deltaT_exch[k] = MedRec.temperature(state_HTF_exch[k]) - MedPB.temperature(state_CO2_exch[k]);
      state_CO2_exch[k].p = p_high;
      state_HTF_exch[k].p = 10 ^ 5;
    end for;
    deltaT_exch[N_exch] = MedRec.temperature(state_HTF_exch[N_exch]) - MedPB.temperature(state_CO2_exch[N_exch]);
    T_HTF_out_des = MedRec.temperature(state_HTF_exch[1]);
    state_CO2_exch[N_exch].p = p_high;
    P_nom = (-turbine.W_turb[1]) - compr.W_comp[1];
    eta_cycle_on = P_nom / (m_des * Q_dis_exch * (N_exch - 1));
//Part I.2. Financial analysis
    C_HTR = if MedPB.temperature(turbine.state_out[1]) >= 550 + 273.15 then 49.45 * sum(HTR.UA_dis) ^ 0.7544 * (1 + 0.02141 * (MedPB.temperature(turbine.state_out[1]) - 550 - 273.15)) else 49.45 * sum(HTR.UA_dis) ^ 0.7544;
    C_turbine = if T_high >= 550 + 273.15 then 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (T_high - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8;
    C_compressor = 1230000 * (compr.W_comp[1] / 10 ^ 6) ^ 0.3392;
    C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
    C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
    C_exchanger = pri_exchanger * m_des * Q_dis_exch * (N_exch - 1) / 1000;
    C_PB = (C_HTR + C_turbine + C_compressor + C_generator + C_cooler + C_exchanger) * 1.05;
//Part II. Off-design points. Mass flow rate is kept constant
    for i in 2:N_points + 1 loop
//    //Inlet of the turbine
      T_high_frac[i - 1] = (0.8 + (1.2 - 0.8) * (i - 1) / (N_points - 1)) * (T_high - 273.15) + 273.15;
      turbine.state_in[i] = MedPB.setState_pTX(p_high, T_high_frac[i - 1]);
//Heat recuperator
      HTR.state_turb[i, N_HTR] = turbine.state_out[i];
      compr.state_in[i] = cooler.state_out;
      HTR.state_comp[i, 1] = compr.state_out[i];
      W_out[i - 1] = (-turbine.W_turb[i]) - compr.W_comp[i];
      eta_cycle_off[i - 1] = W_out[i - 1] / (m_des * (turbine.state_in[i].h - HTR.state_comp[i, N_HTR].h));
      T_CO2_in_HX_table[i - 1] = MedPB.temperature(HTR.state_comp[i, N_HTR]);
    end for;
  equation
    connect(cooler.m_sup, m_sup);
    connect(T_amb, cooler.T_amb);
//Exchanger code
    m_sup = fluid_a.m_flow >= 0.5 * m_des;
    m_HTF_bis = if m_sup then fluid_a.m_flow else m_des;
    for k in 1:N_exch - 1 loop
      m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = m_HTF_bis * (state_HTF_HX[k + 1].h - state_HTF_HX[k].h);
      m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = UA_exch_dis[k] * (deltaT_HX[k] + deltaT_HX[k + 1]) / 2;
      deltaT_HX[k] = MedRec.temperature(state_HTF_HX[k]) - MedPB.temperature(state_CO2_HX[k]);
      state_HTF_HX[k].p = state_HTF_HX[k + 1].p;
      state_CO2_HX[k].p = state_CO2_HX[k + 1].p;
    end for;
    deltaT_HX[N_exch] = MedRec.temperature(state_HTF_HX[N_exch]) - MedPB.temperature(state_CO2_HX[N_exch]);
    state_HTF_HX[N_exch] = if m_sup then MedRec.setState_phX(fluid_a.p, inStream(fluid_a.h_outflow)) else state_HTF_exch[N_exch];
//The TIT and T_CO2 at the inlet of the HX are interpolated from previous points
    T_HTF_in = MedRec.temperature(state_HTF_HX[N_exch]);
    T_HTF_out = MedRec.temperature(state_HTF_HX[1]);
    TIT_CO2 = MedPB.temperature(state_CO2_HX[N_exch]);
    T_CO2_in_HX = Modelica.Math.Vectors.interpolate(T_high_frac, T_CO2_in_HX_table, TIT_CO2);
    state_CO2_HX[1] = MedPB.setState_pTX(p_high, T_CO2_in_HX);
// Power Block calculations
    W_net = if m_sup then max(0, Modelica.Math.Vectors.interpolate(T_high_frac, W_out, TIT_CO2) - cooler.P_cooling) else 0;
//Connectors obligations
    fluid_b.p = fluid_a.p;
    fluid_a.m_flow + fluid_b.m_flow = 0;
    fluid_b.h_outflow = state_HTF_HX[1].h;
    fluid_a.h_outflow = 0;
//shouldn't flow back
    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end interpolWithCompo;
  package ComplicatedPB
  model Turbine "OD model of a turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter SI.Efficiency eta_design = 0.9 "isentropic efficiency of the turbine";
      parameter SI.Efficiency PR = 3 "Pressure ratio";
      parameter SI.AngularVelocity N_shaft = 3358;
      parameter SI.Power P_nom = 100 * 10 ^ 6;
      parameter SI.AbsolutePressure p_high = 250 * 10 ^ 5;
      //parameters for on-design
      parameter Modelica.SIunits.Area A_nozzle(fixed = false);
      parameter SI.Diameter diam_turb(fixed = false);
      parameter SI.Velocity tipSpeed(fixed = false);
      //parameters for off and on-design
      parameter Integer N_points = 10;
      parameter SI.MassFlowRate[N_points + 1] m_flow(each fixed = false, each start = P_nom / 10 ^ 5);
      parameter SI.Velocity[N_points + 1] C_spouting(each fixed = false, each start = 600);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false, each start = 2900) "inlet entropy to the turbine";
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false, each p.start = p_high, each h.start = 1.2 * 10 ^ 6);
      parameter MedPB.ThermodynamicState[N_points + 1] state_isen(each p.fixed = false, each h.fixed = false, each p.start = p_high / PR);
      parameter SI.AbsolutePressure[N_points] p_in(each fixed = false, each start = p_high);
      parameter SI.ThermodynamicTemperature[N_points] T_in(each fixed = false);
      parameter SI.AbsolutePressure[N_points] p_out(each fixed = false, each start = p_high / PR);
      parameter SI.Power[N_points + 1] W_turb(each fixed = false);
      parameter SI.Efficiency[N_points] eta_turb(each fixed = false, each start = eta_design);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
    initial equation
//Part I. On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_isen[1] = MedPB.setState_psX(p_high / PR, s_in[1]);
      state_out[1] = MedPB.setState_phX(p_high / PR, state_in[1].h + eta_design * (state_isen[1].h - state_in[1].h));
      W_turb[1] = m_flow[1] * (state_out[1].h - state_in[1].h);
      C_spouting[1] ^ 2 = 2 * (state_in[1].h - state_isen[1].h);
      m_flow[1] = C_spouting[1] * A_nozzle * MedPB.density(state_out[1]);
      tipSpeed = N_shaft * diam_turb / 2;
      tipSpeed / C_spouting[1] = 0.707;
//Part II. Off-design
      for i in 2:N_points + 1 loop
        state_in[i] = MedPB.setState_pTX(p_in[i - 1], T_in[i - 1]);
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_isen[i] = MedPB.setState_psX(p_out[i - 1], s_in[i]);
        state_out[i] = MedPB.setState_phX(p_out[i - 1], state_in[i].h + eta_turb[i - 1] * (state_isen[i].h - state_in[i].h));
        C_spouting[i] ^ 2 = 2 * (state_in[i].h - state_isen[i].h);
        m_flow[i] = C_spouting[i] * A_nozzle * MedPB.density(state_out[i]);
        eta_turb[i - 1] = eta_design * 2 * (tipSpeed / C_spouting[i]) * sqrt(1 - (tipSpeed / C_spouting[i]) ^ 2);
//eta_turb[i-1]=eta_design;
        W_turb[i] = m_flow[i] * (state_out[i].h - state_in[i].h);
      end for;
      annotation(
        Documentation(info = "<html>
    <p>This turbine's model is based on the phD thesis of J. Dyreby.&nbsp;</p>
  <p>The isentropic efficiency is calculated as a function of the tip speed ration between the tip speed of the rotor and the spouting velocity. It is said to be functionnal for any size.</p>
  <p>The outlet pressure goes beyond the critical pressure for a mass flow too small. The cycle calculation should therefore not be performed below this pressure.</p>
  <p>J. J. Dyreby, &laquo; Modeling the supercritical carbon dioxide Brayton cycle with recompression &raquo;, The University of Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
    </html>"));
      annotation(
        Diagram(graphics = {Text(origin = {-36, -28}, extent = {{18, 80}, {78, 16}}, textString = "TURBINE"), Polygon(origin = {15, 20}, points = {{-35, 44}, {-35, -52}, {35, -68}, {35, 68}, {-35, 44}, {35, 68}, {-35, 44}})}, coordinateSystem(initialScale = 0.1)),
  Icon(graphics = {Text(origin = {-10, 26}, extent = {{-10, 12}, {52, -34}}, textString = "TURBINE"), Ellipse(extent = {{56, 58}, {56, 58}}, endAngle = 360), Polygon(origin = {11, 17}, points = {{-37, 49}, {-37, -51}, {37, -71}, {37, 71}, {-37, 49}})}, coordinateSystem(initialScale = 0.1)));
    end Turbine;

    model Cooler
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      input SI.ThermodynamicTemperature T_amb "Ambiant temperature in Kelvin";
      SI.Power P_cooling "Cooling power necessary to cool down the fluid";
      parameter SI.ThermodynamicTemperature T_amb_des = 40 + 273.15 "Ambiant temperature in Kelvin at design point";
      parameter SI.ThermodynamicTemperature T_low = 55 + 273.15;
      parameter SI.Power P_nom_des = 164000;
      parameter Integer N_cooler = 15;
      //On-design only
      parameter MedPB.ThermodynamicState[N_cooler] state_cooler(each p.fixed = false, each h.fixed = false, each h.start = 4.5 * 10 ^ 5);
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SI.Power P_cool_des(fixed = false) "on-design power necessary to run the fans";
      parameter SI.HeatFlowRate Q_dis_cooler(fixed = false);
      parameter SolarTherm.Types.Conductance[N_cooler - 1] UA_cooler_dis(each fixed = false);
      parameter SolarTherm.Types.Conductance UA_cooler(fixed = false) "Conductance of the cooler in W/K";
      parameter MedPB.ThermodynamicState state_out(p.fixed = false, h.fixed = false);
    initial equation
      for k in 1:N_cooler - 1 loop
        Q_dis_cooler = state_cooler[k + 1].h - state_cooler[k].h;
        m_des * Q_dis_cooler = -UA_cooler_dis[k] * (MedPB.temperature(state_cooler[k + 1]) - T_amb_des + MedPB.temperature(state_cooler[k]) - T_amb_des) / 2;
        if k > 1 then
          state_cooler[k].p = state_cooler[1].p;
        end if;
      end for;
      state_cooler[N_cooler] = state_out;
      UA_cooler = sum(UA_cooler_dis);
      P_cool_des * (T_low - T_amb_des) / (-m_des * Q_dis_cooler * (N_cooler - 1)) = 1.49 * 10 ^ 6 * (35.7 - 30) / (136.6 * 10 ^ 6);
    equation
      P_cooling = P_cool_des * ((T_low - T_amb_des) / (max(T_amb + 5, T_low) - T_amb)) ^ (3 / 0.805);
      annotation(
        Documentation(info = "<html>
    <p>The cooler is thought to be a dry-air cooling device. The outlet temperature of the CO2 is imposed as max(T_low_cycle,T_amb+3). The variation of the ambiant temperature is taken into account in the estimation of the electricity demand for the fans, such as: P_cooling*deltaT/Q_cooler is a constant, deltaT being the average of the temperature of the CO2 and the ambiant, and Q_cooler the energy to withdraw.</p>
    </html>"));
      annotation(
        Icon(graphics = {Rectangle(origin = {2, 1}, extent = {{-58, 65}, {58, -65}}), Text(origin = {0, -1}, extent = {{-40, -15}, {40, 15}}, textString = "COOLER")}),
  Diagram(graphics = {Rectangle(origin = {-4, 7}, extent = {{-64, 67}, {64, -67}}), Text(origin = {5, 14}, extent = {{-41, -12}, {41, 12}}, textString = "COOLER")}));
    end Cooler;

    model HeatRecuperatorDTAve "The heat recuperator is subdivised in N_q segments in order to accurately represent the CO2 properties variation."
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter Integer N_q = 15 "Number of subdivision of the HX";
      parameter Real ratio_m_des = 1 "ratio of m_comp_des/m_turb_des; we suppose m_turb_des=1, and then scale-up";
      parameter Real pinchRecuperator = 5 "pinch of the recuperator. Imposed as a closing equation for on-design";
      parameter SI.Power P_nom = 100 * 10 ^ 6 "on-design power, useful for start values";
      //On-design
      parameter SI.MassFlowRate[N_points + 1] m_flow_turb(each fixed = false);
      parameter SI.MassFlowRate[N_points + 1] m_flow_comp(each fixed = false);
      parameter SolarTherm.Types.Conductance[N_q - 1] UA_dis(each fixed = false, each start = P_nom * 0.4 / N_q);
      parameter SI.HeatFlowRate Q_dis(fixed = false, start = 2 * 10 ^ 4);
      //Off-design
      parameter Integer N_points = 10;
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_turb(each p.fixed = false, each h.fixed = false, each h.start = 7 * 10 ^ 5);
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_comp(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
      parameter SI.TemperatureDifference[N_points + 1, N_q] deltaT(each start = 20, each fixed = false);
    initial equation
//On-design
      for k in 1:N_q - 1 loop
        Q_dis = state_turb[1, k + 1].h - state_turb[1, k].h;
        Q_dis = ratio_m_des * (state_comp[1, k + 1].h - state_comp[1, k].h);
        deltaT[1, k] = MedPB.temperature(state_turb[1, k]) - MedPB.temperature(state_comp[1, k]);
        m_flow_turb[1] * Q_dis = UA_dis[k] * (deltaT[1, k] + deltaT[1, k + 1]) / 2;
        state_comp[1, k + 1].p = state_comp[1, 1].p;
        state_turb[1, k].p = state_turb[1, N_q].p;
      end for;
      deltaT[1, N_q] = MedPB.temperature(state_turb[1, N_q]) - MedPB.temperature(state_comp[1, N_q]);
      min(deltaT[1]) = pinchRecuperator;
//Off-design
      for i in 2:N_points + 1 loop
        for k in 1:N_q - 1 loop
          state_turb[i, k + 1].h - state_turb[i, k].h = ratio_m_des * (state_comp[i, k + 1].h - state_comp[i, k].h);
          deltaT[i, k] = MedPB.temperature(state_turb[i, k]) - MedPB.temperature(state_comp[i, k]);
          m_flow_turb[i] * (state_turb[i, k + 1].h - state_turb[i, k].h) = UA_dis[k] * (deltaT[i, k] + deltaT[i, k + 1]) / 2;
          state_turb[i, k].p = state_turb[i, N_q].p;
          state_comp[i, k + 1].p = state_comp[i, 1].p;
        end for;
        deltaT[i, N_q] = MedPB.temperature(state_turb[i, N_q]) - MedPB.temperature(state_comp[i, N_q]);
      end for;
      annotation(
        Diagram(graphics = {Rectangle(origin = {1, 7}, extent = {{-61, 31}, {61, -31}}), Text(origin = {5, 1}, extent = {{-53, -17}, {53, 17}}, textString = "RECUPERATOR")}),
  Icon(graphics = {Rectangle(origin = {-3, -9}, extent = {{-65, 33}, {65, -33}}), Text(origin = {-2, -5}, extent = {{-46, -15}, {46, 15}}, textString = "RECUPERATOR")}));
      annotation(
        Documentation(info = "<html>
    <p>This heat recuperator is a counter-flow HX. Closure equations are based on the equality of m_flow*delta_H for both sides and m_flow*delta_H= UA_i*DTAve_i, DTAve being the average of the temperature difference between the inlet and the outlet of the sub-HX.</p>
  <p>The UA_i must be given as parameters from the on-design analysis.&nbsp;</p>
    
    </html>"));
    end HeatRecuperatorDTAve;

    model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      parameter SI.Efficiency eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
      parameter SI.AngularVelocity N_design = 3358 "Design rotationnal speed in rad/s";
      parameter Real PR = 3 "pressure ratio chosen";
      //On-design calculations
      parameter SI.AbsolutePressure p_high = CV.from_bar(250);
      parameter SI.MassFlowRate[N_points + 1] m_flow(each fixed = false);
      parameter SI.Diameter diam_rotor(fixed = false) "on-design diameter of the rotor";
      parameter SI.Velocity tipSpeed(fixed = false) " tip speed of the rotor";
      parameter Real psi_des(fixed = false) "on-design adimensionned head";
      parameter Real phi_opt = 0.0297035 "optimal adimensionned mass flow";
      //Off design calculations
      parameter Integer N_points = 10;
      parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false, each h.start = 450000, each p.start = p_high / PR);
      parameter MedPB.ThermodynamicState[N_points + 1] state_isen(each p.fixed = false, each h.fixed = false, each p.start = p_high, each h.start = 500000);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 550000, each p.start = p_high);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false, each start = 1500);
      parameter Real[N_points] psi(each fixed = false, each start = psi_des) "off-design adimensionned head";
      parameter Real[N_points] phi(each fixed = false, each start = phi_opt) "off-design adimensionned mass flow";
      parameter SI.AbsolutePressure[N_points] p_out(each fixed = false, each start = p_high);
      parameter SI.Efficiency[N_points] eta_comp(each fixed = false, each start = eta_design);
      parameter SI.Power[N_points + 1] W_comp(each fixed = false);
    initial equation
//On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_isen[1] = MedPB.setState_psX(state_in[1].p * PR, s_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p * PR, state_in[1].h + (state_isen[1].h - state_in[1].h) / eta_design);
      W_comp[1] = m_flow[1] * (state_out[1].h - state_in[1].h);
      2 * m_flow[1] = phi_opt * MedPB.density(state_in[1]) * N_design * diam_rotor ^ 3;
      tipSpeed = N_design * diam_rotor / 2;
      psi_des = (state_isen[1].h - state_in[1].h) / tipSpeed ^ 2;
//Off-design
      for i in 2:N_points + 1 loop
        state_in[i] = MedPB.setState_pTX(p_high / PR, T_low);
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_isen[i] = MedPB.setState_psX(p_out[i - 1], s_in[i]);
        state_out[i] = MedPB.setState_phX(p_out[i - 1], state_in[i].h + (state_isen[i].h - state_in[i].h) / eta_comp[i - 1]);
        W_comp[i] = m_flow[i] * (state_out[i].h - state_in[i].h);
        phi[i - 1] = m_flow[i] / (MedPB.density(state_in[i]) * tipSpeed * diam_rotor ^ 2);
        psi[i - 1] = (state_isen[i].h - state_in[i].h) / tipSpeed ^ 2;
        psi[i - 1] = (0.04049 + 54.7 * phi[i - 1] - 2505 * phi[i - 1] ^ 2 + 53224 * phi[i - 1] ^ 3 - 498626 * phi[i - 1] ^ 4) * psi_des / 0.46181921979961293;
        eta_comp[i - 1] = eta_design / 0.677837 * ((-0.7069) + 168.6 * phi[i - 1] - 8089 * phi[i - 1] ^ 2 + 182725 * phi[i - 1] ^ 3 - 1.638 * 10 ^ 6 * phi[i - 1] ^ 4);
//eta_comp[i-1]=eta_design;
      end for;
    equation

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

    model simpleRecupPB "This model calculates a few off-design points as initial equation and after interpolate over them"
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      import FI = SolarTherm.Models.Analysis.Finances;
      extends SolarTherm.Media.CO2.PropCO2;
      extends Icons.PowerBlock;
      input SI.ThermodynamicTemperature T_amb;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
      //Objective: having a few points in order to interpolate the HTF outlet temperature and the TIT. Should be feasible with for loops in the initial equation.
      parameter SI.Power P_nom = 100 * 10 ^ 6;
      parameter SI.AbsolutePressure p_high = CV.from_bar(200);
      parameter SI.AbsolutePressure p_low = CV.from_bar(80);
      parameter Real PR = p_high / p_low;
      parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_HTF_in_des = CV.from_degC(800) "inlet temperature of the HTF";
      parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
      parameter SI.ThermodynamicTemperature T_amb_des = CV.from_degC(20) "ambient temperature at design";
      parameter SI.MassFlowRate m_des(fixed = false, start = P_nom / 10 ^ 5);
      parameter SI.TemperatureDifference pinchHTR = 5;
      parameter SI.TemperatureDifference pinchExch = T_HTF_in_des - T_high;
      parameter SI.Efficiency eta_comp = 0.89;
      parameter SI.Efficiency eta_turb = 0.93;
      parameter Integer N_exch = 5;
      parameter Integer N_HTR = 15;
      parameter Integer N_cooler = 15;
      //Variables for real time calculation
      //Exchanger variables
      Boolean m_sup "indicates if mass flow from tank is superior enough to calculate";
      SI.MassFlowRate m_HTF_bis "mass flow rate used for switching off the PB";
      SI.ThermodynamicTemperature T_HTF_in "Inlet temperature of HTF";
      SI.ThermodynamicTemperature T_HTF_out "Outlet temperature of HTF";
      SI.ThermodynamicTemperature TIT_CO2 "Turbine inlet Temperature of CO2";
      SI.ThermodynamicTemperature T_CO2_in_HX "CO2 temperature at the inlet of the exchanger";
      MedPB.ThermodynamicState[N_exch] state_CO2_HX(each h.start = 10 ^ 6);
      MedRec.ThermodynamicState[N_exch] state_HTF_HX(each h.start = 6 * 10 ^ 5);
      SI.TemperatureDifference[N_exch] deltaT_HX(each start = 20);
      SI.MassFlowRate m_CO2(start = m_des);
      //PowerBlock variables
      Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
        Placement(transformation(extent = {{-54, 22}, {-34, 42}}), iconTransformation(extent = {{-48, 30}, {-40, 38}})));
      Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
        Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
      Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
        Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
      parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
      //Exchanger parameters. Calculated only at on-design
      parameter SolarTherm.Types.Conductance[N_exch - 1] UA_exch_dis(each fixed = false, each start = P_nom * 0.4 / N_exch);
      parameter SI.HeatFlowRate Q_dis_exch(fixed = false);
      parameter MedPB.ThermodynamicState[N_exch] state_CO2_exch(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
      parameter MedRec.ThermodynamicState[N_exch] state_HTF_exch(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
      parameter SI.TemperatureDifference[N_exch] deltaT_exch(each start = 20, each fixed = false);
      parameter SI.ThermodynamicTemperature T_HTF_out_des(fixed = false);
      //Off-design parameters for calculation
      parameter Integer N_points = 30;
      parameter SI.ThermodynamicTemperature[N_points] T_high_frac(each fixed = false) "TIT points for off-design";
      parameter SI.ThermodynamicTemperature[N_points] T_CO2_in_HX_table(each fixed = false) "Temperature at the inlet of the exchanger, for a TIT given. Used for interpolation";
      parameter SI.MassFlowRate[N_points] m_flow(each fixed = false, each start = m_des) "mass flow rate in the cycle for a given off-design point";
      parameter SI.Power[N_points] W_out(each fixed = false) "Outlet power at off-design";
      parameter SI.Efficiency eta_cycle_on(fixed = false, start = 0.5) "efficiency of the cycle";
      //parameter SI.Efficiency[N_points] eta_cycle_off(each fixed = false, each start = 0.5) "efficiency of the cycle";
      //Financial analysis
      parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
      parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
      parameter FI.Money C_compressor(fixed = false) "cost of the main compressor";
      parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
      parameter FI.Money C_generator(fixed = false) "cost of the generator";
      parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
      parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.ComplicatedPB.Turbine turbine(eta_design = eta_turb, PR = PR, N_points = N_points, p_high = p_high) annotation(
        Placement(visible = true, transformation(origin = {38, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.ComplicatedPB.Cooler cooler annotation(
        Placement(visible = true, transformation(origin = {-52, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.ComplicatedPB.HeatRecuperatorDTAve HTR(N_points = N_points, N_q = N_HTR, pinchRecuperator = pinchHTR, ratio_m_des = 1, P_nom = P_nom) annotation(
        Placement(visible = true, transformation(origin = {-2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.ComplicatedPB.CompressorOnShaft compr(eta_design = eta_comp, PR = PR, N_points = N_points, p_high = p_high) annotation(
        Placement(visible = true, transformation(origin = {-38, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
//Part I. On-design=first point
//Inlet of the turbine
      turbine.state_in[1] = MedPB.setState_pTX(p_high, T_high);
      m_des = turbine.m_flow[1];
      turbine.m_flow[1] = HTR.m_flow_turb[1];
//High temperature recuperator
      HTR.state_turb[1, N_HTR] = turbine.state_out[1];
      cooler.state_cooler[1] = HTR.state_turb[1, 1];
      cooler.m_des = HTR.m_flow_turb[1];
//Cooler. UA (and therefore this part) are only for economical correlation
      cooler.state_out = MedPB.setState_pTX(p_low, T_low);
//Main compressor. Inlet supposed always the same for off-design
      compr.state_in[1] = cooler.state_out;
      compr.m_flow[1] = cooler.m_des;
      HTR.state_comp[1, 1] = compr.state_out[1];
      HTR.m_flow_comp[1] = compr.m_flow[1];
//Exchanger. Done only on-design. After, only the TIT is varied.
      state_CO2_exch[N_exch].h = turbine.state_in[1].h;
      state_CO2_exch[1].h = HTR.state_comp[1, N_HTR].h;
      state_HTF_exch[N_exch] = MedRec.setState_pTX(10 ^ 5, T_HTF_in_des);
      for k in 1:N_exch - 1 loop
        Q_dis_exch = state_CO2_exch[k + 1].h - state_CO2_exch[k].h;
        Q_dis_exch = state_HTF_exch[k + 1].h - state_HTF_exch[k].h;
        m_des * Q_dis_exch = UA_exch_dis[k] * (deltaT_exch[k] + deltaT_exch[k + 1]) / 2;
        deltaT_exch[k] = MedRec.temperature(state_HTF_exch[k]) - MedPB.temperature(state_CO2_exch[k]);
        state_CO2_exch[k].p = p_high;
        state_HTF_exch[k].p = 10 ^ 5;
      end for;
      deltaT_exch[N_exch] = MedRec.temperature(state_HTF_exch[N_exch]) - MedPB.temperature(state_CO2_exch[N_exch]);
      T_HTF_out_des = MedRec.temperature(state_HTF_exch[1]);
      state_CO2_exch[N_exch].p = p_high;
      P_nom = (-turbine.W_turb[1]) - compr.W_comp[1];
      eta_cycle_on = P_nom / (m_des * Q_dis_exch * (N_exch - 1));
//Part I.2. Financial analysis
      C_HTR = if MedPB.temperature(turbine.state_out[1]) >= 550 + 273.15 then 49.45 * sum(HTR.UA_dis) ^ 0.7544 * (1 + 0.02141 * (MedPB.temperature(turbine.state_out[1]) - 550 - 273.15)) else 49.45 * sum(HTR.UA_dis) ^ 0.7544;
      C_turbine = if T_high >= 550 + 273.15 then 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (T_high - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8;
      C_compressor = 1230000 * (compr.W_comp[1] / 10 ^ 6) ^ 0.3392;
      C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
      C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
      C_exchanger = pri_exchanger * m_des * Q_dis_exch * (N_exch - 1) / 1000;
      C_PB = (C_HTR + C_turbine + C_compressor + C_generator + C_cooler + C_exchanger) * 1.05;
//Part II. Off-design points. Mass flow rate is kept constant
      for i in 2:N_points + 1 loop
//Inlet of the turbine
        T_high_frac[i - 1] = (0.8 + (1.2 - 0.8) * (i - 2) / (N_points - 1)) * (T_high - 273.15) + 273.15;
        turbine.T_in[i - 1] = T_high_frac[i - 1];
        turbine.p_in[i - 1] = p_high;
        turbine.p_out[i - 1] = turbine.p_in[i - 1] / PR;
        m_flow[i - 1] = turbine.m_flow[i];
        turbine.m_flow[i] = HTR.m_flow_turb[i];
//Heat recuperator
        HTR.state_turb[i, N_HTR] = turbine.state_out[i];
        HTR.state_comp[i, 1] = compr.state_out[i];
        HTR.m_flow_turb[i] = compr.m_flow[i];
        compr.m_flow[i] = HTR.m_flow_comp[i];
        W_out[i - 1] = (-turbine.W_turb[i]) - compr.W_comp[i];
//eta_cycle_off[i - 1] = W_out[i - 1] / (m_flow[i-1] * (turbine.state_in[i].h - HTR.state_comp[i, N_HTR].h));
        T_CO2_in_HX_table[i - 1] = MedPB.temperature(HTR.state_comp[i, N_HTR]);
      end for;
    equation
      connect(T_amb, cooler.T_amb);
//Exchanger code
      m_sup = fluid_a.m_flow >= 0.5 * m_des;
      m_HTF_bis = if m_sup then fluid_a.m_flow else m_des;
      m_CO2 = if m_sup then Modelica.Math.Vectors.interpolate(T_high_frac, m_flow, TIT_CO2) else m_des;
      for k in 1:N_exch - 1 loop
        m_CO2 * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = m_HTF_bis * (state_HTF_HX[k + 1].h - state_HTF_HX[k].h);
        m_CO2 * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = UA_exch_dis[k] * (deltaT_HX[k] + deltaT_HX[k + 1]) / 2;
        deltaT_HX[k] = MedRec.temperature(state_HTF_HX[k]) - MedPB.temperature(state_CO2_HX[k]);
        state_HTF_HX[k].p = state_HTF_HX[k + 1].p;
        state_CO2_HX[k].p = state_CO2_HX[k + 1].p;
      end for;
      deltaT_HX[N_exch] = MedRec.temperature(state_HTF_HX[N_exch]) - MedPB.temperature(state_CO2_HX[N_exch]);
      state_HTF_HX[N_exch] = if m_sup then MedRec.setState_phX(fluid_a.p, inStream(fluid_a.h_outflow)) else state_HTF_exch[N_exch];
//The TIT and T_CO2 at the inlet of the HX are interpolated from previous points
      T_HTF_in = MedRec.temperature(state_HTF_HX[N_exch]);
      T_HTF_out = MedRec.temperature(state_HTF_HX[1]);
      TIT_CO2 = MedPB.temperature(state_CO2_HX[N_exch]);
      T_CO2_in_HX = Modelica.Math.Vectors.interpolate(T_high_frac, T_CO2_in_HX_table, TIT_CO2);
      state_CO2_HX[1] = MedPB.setState_pTX(p_high, T_CO2_in_HX);
// Power Block calculations
      W_net = if m_sup then max(0, Modelica.Math.Vectors.interpolate(T_high_frac, W_out, TIT_CO2) - cooler.P_cooling) else 0;
//Connectors obligations
      fluid_b.p = fluid_a.p;
      fluid_a.m_flow + fluid_b.m_flow = 0;
      fluid_b.h_outflow = state_HTF_HX[1].h;
      fluid_a.h_outflow = 0;
//shouldn't flow back
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end simpleRecupPB;

    model test
      extends SolarTherm.Media.CO2.PropCO2;
      import CV = Modelica.SIunits.Conversions;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
      SourceFlow src(T_out = 800 + 273.15, p_out = 10 ^ 5, m_flow = PB.m_des, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
        Placement(visible = true, transformation(origin = {-28, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {-26, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.ComplicatedPB.interpolWithCompo PB annotation(
        Placement(visible = true, transformation(origin = {31, 15}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    initial equation

    equation
      connect(src.port_b, PB.fluid_a) annotation(
        Line(points = {{-20, 54}, {18, 54}, {18, 26}, {18, 26}, {18, 24}}, color = {0, 127, 255}));
      connect(sink.port_a, PB.fluid_b) annotation(
        Line(points = {{-18, 12}, {-16, 12}, {-16, 2}, {14, 2}, {14, 2}}, color = {0, 127, 255}));
      PB.T_amb = 273.15;
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end test;
  end ComplicatedPB;

  package Simpler
  model Turbine "OD model of a turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter SI.Efficiency eta_design = 0.9 "isentropic efficiency of the turbine";
      parameter SI.Efficiency PR = 3 "Pressure ratio";
      parameter SI.AngularVelocity N_shaft = 3358;
      parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "TIT at design point";
      //parameters for off and on-design
      parameter Integer N_points = 10;
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false, each start = 2900) "inlet entropy to the turbine";
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false);
      parameter SI.Power[N_points + 1] W_turb(each fixed = false);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false);
    initial equation
//Part I. On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p / PR, state_in[1].h + eta_design * (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[1].p / PR, s_in[1])) - state_in[1].h));
      W_turb[1] = m_des * (state_out[1].h - state_in[1].h);
      for i in 2:N_points + 1 loop
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_out[i] = MedPB.setState_phX(state_in[i].p / PR, state_in[i].h + eta_design * (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[i].p / PR, s_in[i])) - state_in[i].h));
        W_turb[i] = m_des * (state_out[i].h - state_in[i].h);
      end for;
      annotation(
        Documentation(info = "<html>
  		<p>This turbine's model is based on the phD thesis of J. Dyreby.&nbsp;</p>
  <p>The isentropic efficiency is calculated as a function of the tip speed ration between the tip speed of the rotor and the spouting velocity. It is said to be functionnal for any size.</p>
  <p>The outlet pressure goes beyond the critical pressure for a mass flow too small. The cycle calculation should therefore not be performed below this pressure.</p>
  <p>J. J. Dyreby, &laquo; Modeling the supercritical carbon dioxide Brayton cycle with recompression &raquo;, The University of Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
  		</html>"));
      annotation(
        Diagram(graphics = {Text(origin = {-36, -28}, extent = {{18, 80}, {78, 16}}, textString = "TURBINE"), Polygon(origin = {15, 20}, points = {{-35, 44}, {-35, -52}, {35, -68}, {35, 68}, {-35, 44}, {35, 68}, {-35, 44}})}, coordinateSystem(initialScale = 0.1)),
  Icon(graphics = {Text(origin = {-10, 26}, extent = {{-10, 12}, {52, -34}}, textString = "TURBINE"), Ellipse(extent = {{56, 58}, {56, 58}}, endAngle = 360), Polygon(origin = {11, 17}, points = {{-37, 49}, {-37, -51}, {37, -71}, {37, 71}, {-37, 49}})}, coordinateSystem(initialScale = 0.1)));
    end Turbine;

    model Exchanger
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
      import SI = Modelica.SIunits;
      parameter SI.ThermodynamicTemperature T_out_CO2_des = 715 + 273.15;
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
      SI.SpecificEnthalpy[N_exch] h_CO2(start = {990000 + i / N_exch * 200000 for i in 1:N_exch}), h_HTF(start = {600000 + i / N_exch * 200000 for i in 1:N_exch});
      Real[N_exch] deltaT "Temperature difference in the heat exchangers";
      SI.HeatFlowRate Q_HX;
      SI.ThermodynamicTemperature T_CO2_out, T_HTF_out;
      //Real deltaT_lm;
      Real deltaTAve;
      SI.MassFlowRate m_HTF_bis(start = P_nom_des / 10 ^ 5);
      parameter SI.HeatFlowRate Q_HX_des(fixed = false);
      parameter SI.MassFlowRate m_CO2_des(fixed = false), m_HTF_des(fixed = false);
      parameter SolarTherm.Types.Conductance UA_HX(fixed = false) "on-design conductance of the overall exchanger";
      parameter SolarTherm.Types.Conductance[N_exch - 1] UA_HX_dis(each fixed = false) "on-design conductance of the exchanger";
      parameter SI.SpecificEnthalpy h_in_HTF_des(fixed = false, start = 855000), h_out_HTF_des(fixed = false), h_in_CO2_des(fixed = false, start = 900000), h_out_CO2_des(fixed = false, start = 1.2 * 10 ^ 6);
      parameter Real[N_exch] deltaT_des(each fixed = false, each start = 75);
      parameter SI.AbsolutePressure p_in_CO2_des(fixed = false), p_out_CO2_des(fixed = false);
      parameter SI.AbsolutePressure p_in_HTF_des(fixed = false), p_out_HTF_des(fixed = false);
      parameter SI.ThermodynamicTemperature[N_exch] T_CO2_des(each fixed = false, start = {600 + 273.15 + 120 * (i / N_exch) for i in 1:N_exch}), T_HTF_des(each fixed = false, start = {650 + 273.15 + 120 * (i / N_exch) for i in 1:N_exch});
      parameter MedPB.ThermodynamicState[N_exch] state_CO2_des(each p.fixed = false, each h.fixed = false, each p.start = 250 * 10 ^ 5, each h.start = 10 ^ 6), state_HTF_des(each p.fixed = false, each h.fixed = false, each p.start = 10 ^ 5, each h.start = 855004);
    initial equation
      for i in 1:N_exch loop
        deltaT_des[i] = MedRec.temperature(state_HTF_des[i]) - MedPB.temperature(state_CO2_des[i]);
        state_CO2_des[i] = MedPB.setState_pTX(p_in_CO2_des, T_CO2_des[i]);
        state_HTF_des[i] = MedRec.setState_pTX(p_in_HTF_des, T_HTF_des[i]);
      end for;
      T_CO2_des[N_exch] = T_out_CO2_des;
      for i in 1:N_exch - 1 loop
        Q_HX_des = ratio_m_des * (state_CO2_des[i + 1].h - state_CO2_des[i].h);
        Q_HX_des = state_HTF_des[i + 1].h - state_HTF_des[i].h;
        m_HTF_des * Q_HX_des = UA_HX_dis[i] * (deltaT_des[i] + deltaT_des[i + 1]) / 2;
      end for;
      UA_HX = sum(UA_HX_dis);
      p_in_CO2_des = p_out_CO2_des;
      p_in_HTF_des = p_out_HTF_des;
      h_in_HTF_des = MedRec.specificEnthalpy(state_HTF_des[N_exch]);
      h_out_HTF_des = MedRec.specificEnthalpy(state_HTF_des[1]);
      h_in_CO2_des = state_CO2_des[1].h;
      h_out_CO2_des = state_CO2_des[N_exch].h;
      m_CO2_des = ratio_m_des * m_HTF_des;
    equation
      for i in 1:N_exch loop
        deltaT[i] = if m_sup then MedRec.temperature(state_HTF[i]) - MedPB.temperature(state_CO2[i]) else deltaT_des[i];
        state_CO2[i] = MedPB.setState_phX(CO2_port_a.p, h_CO2[i]);
        state_HTF[i] = MedRec.setState_phX(HTF_port_a.p, h_HTF[i]);
      end for;
      T_CO2_out = MedPB.temperature(state_CO2[N_exch]);
      T_HTF_out = MedRec.temperature(state_HTF[1]);
//deltaT_lm = if deltaT[2] * deltaT[1] < 0 then (abs(deltaT[1]) ^ (1 / 3) * sign(deltaT[1]) / 2 + abs(deltaT[2]) ^ (1 / 3) * sign(deltaT[2]) / 2) ^ 3 else (deltaT[1] - deltaT[2]) / (Modelica.Math.log(deltaT[1] / deltaT[2]) + 0.0001);
      deltaTAve = (deltaT[1] + deltaT[N_exch]) / 2;
      h_CO2[N_exch] = CO2_port_b.h_outflow;
      h_HTF[N_exch] = if m_sup then inStream(HTF_port_a.h_outflow) else h_in_HTF_des;
      h_CO2[1] = inStream(CO2_port_a.h_outflow);
      HTF_port_b.h_outflow = if m_sup then h_HTF[1] else inStream(HTF_port_a.h_outflow);
      m_HTF_bis = if m_sup then HTF_port_a.m_flow else m_HTF_des;
      Q_HX = CO2_port_a.m_flow * (h_CO2[N_exch] - h_CO2[1]);
      for i in 1:N_exch - 1 loop
        m_HTF_bis * (h_HTF[i + 1] - h_HTF[i]) = CO2_port_a.m_flow * (h_CO2[i + 1] - h_CO2[i]);
        CO2_port_a.m_flow * (h_CO2[i + 1] - h_CO2[i]) = UA_HX_dis[i] * (1 / 2 * abs(m_HTF_bis / m_HTF_des + CO2_port_a.m_flow / m_CO2_des)) ^ 0.8 * (deltaT[i] + deltaT[i + 1]) / 2;
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
        Documentation(info = "<html>
  		<p>The exchanger is a heat exchanger between the HTF and the CO2. It is a counterflow HX, based on a TLMD. The conductance UA has to be specified from the on-design.</p>
  <p>The conductance in off-design varies as UA_Off=UA_on*(m_flow/m_design)^0.8.&nbsp;<span >The average between the two mass flows is taken.</span></p>
  <p>A.T. Louis et T. Neises, analysis and optimization for Off-design performance of the recompression s-CO2 cycles for high temperature CSP applications, in The 5th International Symposium-Supercritical CO2 Power Cycles, 2016</p>
  <p>&nbsp;</p>
  		</html>"));
      annotation(
        Diagram(graphics = {Rectangle(origin = {1, 4}, extent = {{-57, 40}, {57, -40}}), Text(origin = {-1, 8}, extent = {{-47, 16}, {47, -16}}, textString = "Exchanger")}),
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end Exchanger;

    model Cooler
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      input Boolean m_sup;
      input SI.ThermodynamicTemperature T_amb "Ambiant temperature in Kelvin";
      SI.Power P_cooling "Cooling power necessary to cool down the fluid";
      parameter SI.ThermodynamicTemperature T_amb_des = 40 + 273.15 "Ambiant temperature in Kelvin at design point";
      parameter SI.ThermodynamicTemperature T_low = 55 + 273.15;
      parameter SI.Power P_nom_des = 164000;
      parameter Integer N_cooler = 15;
      //On-design only
      parameter MedPB.ThermodynamicState[N_cooler] state_cooler(each p.fixed = false, each h.fixed = false, each h.start = 4.5 * 10 ^ 5);
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SI.Power P_cool_des(fixed = false) "on-design power necessary to run the fans";
      parameter SI.HeatFlowRate Q_dis_cooler(fixed = false);
      parameter SolarTherm.Types.Conductance[N_cooler - 1] UA_cooler_dis(each fixed = false);
      parameter SolarTherm.Types.Conductance UA_cooler(fixed = false) "Conductance of the cooler in W/K";
      parameter MedPB.ThermodynamicState state_out(p.fixed = false, h.fixed = false);
    initial equation
      for k in 1:N_cooler - 1 loop
        Q_dis_cooler = state_cooler[k + 1].h - state_cooler[k].h;
        m_des * Q_dis_cooler = -UA_cooler_dis[k] * (MedPB.temperature(state_cooler[k + 1]) - T_amb_des + MedPB.temperature(state_cooler[k]) - T_amb_des) / 2;
        if k > 1 then
          state_cooler[k].p = state_cooler[1].p;
        end if;
      end for;
      state_cooler[N_cooler] = state_out;
      UA_cooler = sum(UA_cooler_dis);
      P_cool_des * (T_low - T_amb_des) / (-m_des * Q_dis_cooler * (N_cooler - 1)) = 1.49 * 10 ^ 6 * (35.7 - 30) / (136.6 * 10 ^ 6);
    equation
      P_cooling = P_cool_des * ((T_low - T_amb_des) / (max(T_amb + 5, T_low) - T_amb)) ^ (3 / 0.805);
      annotation(
        Documentation(info = "<html>
  		<p>The cooler is thought to be a dry-air cooling device. The outlet temperature of the CO2 is imposed as max(T_low_cycle,T_amb+3). The variation of the ambiant temperature is taken into account in the estimation of the electricity demand for the fans, such as: P_cooling*deltaT/Q_cooler is a constant, deltaT being the average of the temperature of the CO2 and the ambiant, and Q_cooler the energy to withdraw.</p>
  		</html>"));
      annotation(
        Icon(graphics = {Rectangle(origin = {2, 1}, extent = {{-58, 65}, {58, -65}}), Text(origin = {0, -1}, extent = {{-40, -15}, {40, 15}}, textString = "COOLER")}),
  Diagram(graphics = {Rectangle(origin = {-4, 7}, extent = {{-64, 67}, {64, -67}}), Text(origin = {5, 14}, extent = {{-41, -12}, {41, 12}}, textString = "COOLER")}));
    end Cooler;

    model HeatRecuperatorDTAve "The heat recuperator is subdivised in N_q segments in order to accurately represent the CO2 properties variation."
      extends SolarTherm.Media.CO2.PropCO2;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      import SI = Modelica.SIunits;
      parameter Integer N_q = 15 "Number of subdivision of the HX";
      parameter Real ratio_m_des = 1 "ratio of m_comp_des/m_turb_des; we suppose m_turb_des=1, and then scale-up";
      parameter Real pinchRecuperator = 5 "pinch of the recuperator. Imposed as a closing equation for on-design";
      parameter SI.Power P_nom = 100 * 10 ^ 6 "on-design power, useful for start values";
      //On-design
      parameter SI.MassFlowRate m_des(fixed = false);
      parameter SolarTherm.Types.Conductance[N_q - 1] UA_dis(each fixed = false, each start = P_nom * 0.4 / N_q);
      parameter SI.HeatFlowRate Q_dis(fixed = false, start = 2 * 10 ^ 4);
      //Off-design
      parameter Integer N_points = 10;
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_turb(each p.fixed = false, each h.fixed = false, each h.start = 7 * 10 ^ 5);
      parameter MedPB.ThermodynamicState[N_points + 1, N_q] state_comp(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
      parameter SI.TemperatureDifference[N_points + 1, N_q] deltaT(each start = 20, each fixed = false);
    initial equation
//On-design
      for k in 1:N_q - 1 loop
        Q_dis = state_turb[1, k + 1].h - state_turb[1, k].h;
        Q_dis = ratio_m_des * (state_comp[1, k + 1].h - state_comp[1, k].h);
        deltaT[1, k] = MedPB.temperature(state_turb[1, k]) - MedPB.temperature(state_comp[1, k]);
        m_des * Q_dis = UA_dis[k] * (deltaT[1, k] + deltaT[1, k + 1]) / 2;
        state_comp[1, k + 1].p = state_comp[1, 1].p;
        state_turb[1, k].p = state_turb[1, N_q].p;
      end for;
      deltaT[1, N_q] = MedPB.temperature(state_turb[1, N_q]) - MedPB.temperature(state_comp[1, N_q]);
      min(deltaT[1]) = pinchRecuperator;
//Off-design
      for i in 2:N_points + 1 loop
        for k in 1:N_q - 1 loop
          state_turb[i, k + 1].h - state_turb[i, k].h = ratio_m_des * (state_comp[i, k + 1].h - state_comp[i, k].h);
          deltaT[i, k] = MedPB.temperature(state_turb[i, k]) - MedPB.temperature(state_comp[i, k]);
          m_des * (state_turb[i, k + 1].h - state_turb[i, k].h) = UA_dis[k] * (deltaT[i, k] + deltaT[i, k + 1]) / 2;
          state_turb[i, k].p = state_turb[i, N_q].p;
          state_comp[i, k + 1].p = state_comp[i, 1].p;
        end for;
        deltaT[i, N_q] = MedPB.temperature(state_turb[i, N_q]) - MedPB.temperature(state_comp[i, N_q]);
      end for;
      annotation(
        Diagram(graphics = {Rectangle(origin = {1, 7}, extent = {{-61, 31}, {61, -31}}), Text(origin = {5, 1}, extent = {{-53, -17}, {53, 17}}, textString = "RECUPERATOR")}),
  Icon(graphics = {Rectangle(origin = {-3, -9}, extent = {{-65, 33}, {65, -33}}), Text(origin = {-2, -5}, extent = {{-46, -15}, {46, 15}}, textString = "RECUPERATOR")}));
      annotation(
        Documentation(info = "<html>
  		<p>This heat recuperator is a counter-flow HX. Closure equations are based on the equality of m_flow*delta_H for both sides and m_flow*delta_H= UA_i*DTAve_i, DTAve being the average of the temperature difference between the inlet and the outlet of the sub-HX.</p>
  <p>The UA_i must be given as parameters from the on-design analysis.&nbsp;</p>
  		
  		</html>"));
    end HeatRecuperatorDTAve;

    model CompressorOnShaft "0D model of a compressor on the same shaft as the turbine"
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      parameter SI.Efficiency eta_design = 0.89 "Maximal isentropic efficiency of the compressor";
      parameter SI.AngularVelocity N_design = 40000 * 0.104 "Design rotationnal speed in rad/s";
      parameter Real PR = 3 "pressure ratio chosen";
      //On-design calculations
      parameter SI.AbsolutePressure p_high = CV.from_bar(250);
      parameter SI.MassFlowRate m_des(fixed = false);
      //Off design calculations
      parameter Integer N_points = 10;
      parameter MedPB.ThermodynamicState[N_points + 1] state_in(each p.fixed = false, each h.fixed = false, each h.start = 550000);
      parameter MedPB.ThermodynamicState[N_points + 1] state_out(each p.fixed = false, each h.fixed = false, each h.start = 550000);
      parameter SI.SpecificEntropy[N_points + 1] s_in(each fixed = false);
      parameter SI.Power[N_points + 1] W_comp(each fixed = false);
    initial equation
//On-design
      s_in[1] = MedPB.specificEntropy(state_in[1]);
      state_out[1] = MedPB.setState_phX(state_in[1].p * PR, state_in[1].h + (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[1].p * PR, s_in[1])) - state_in[1].h) / eta_design);
      W_comp[1] = m_des * (state_out[1].h - state_in[1].h);
//Off-design
      for i in 2:N_points + 1 loop
        s_in[i] = MedPB.specificEntropy(state_in[i]);
        state_out[i] = MedPB.setState_phX(state_in[i].p * PR, state_in[i].h + (MedPB.specificEnthalpy(MedPB.setState_psX(state_in[i].p * PR, s_in[i])) - state_in[i].h) / eta_design);
        W_comp[i] = m_des * (state_out[i].h - state_in[i].h);
      end for;
    equation

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

    model recompPB
      extends SolarTherm.Media.CO2.PropCO2;
      import SI = Modelica.SIunits;
      import FI = SolarTherm.Models.Analysis.Finances;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
      //input SolarTherm.Interfaces.Connectors.WeatherBus wbus;
      extends Icons.PowerBlock;
      //  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
      //    Placement(visible = true,transformation(extent = {{-54, 22}, {-34, 42}}, rotation = 0), iconTransformation(extent = {{-48, 28}, {-40, 36}}, rotation = 0)));
      //  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
      //    Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
      Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
        Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{48, -12}, {58, -2}}, rotation = 0)));
      // PB parameters
      parameter Boolean external_parasities = false "= true enable parasities as an input";
      parameter Real nu_min = 0.25 "Minimum turbine operation";
      Modelica.Blocks.Interfaces.RealInput parasities if external_parasities annotation(
        Placement(transformation(extent = {{-12, -12}, {12, 12}}, rotation = -90, origin = {1.77636e-015, 80}), iconTransformation(extent = {{-6, -6}, {6, 6}}, rotation = -90, origin = {20, 60})));
      input SI.ThermodynamicTemperature T_amb;
      //Cycle parameters
      parameter SI.AbsolutePressure p_high = 200 * 10 ^ 5 "high pressure of the cycle";
      parameter SI.ThermodynamicTemperature T_high = 715 + 273.15 "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_amb_des = 30 + 273.15 "ambiant temperature";
      parameter Real PR = 2.5 "Pressure ratio";
      parameter SI.Power P_gro = 100 * 10 ^ 6 "first guess of power outlet";
      parameter SI.Power P_nom(fixed = false) "Electrical power at design point";
      parameter SI.MassFlowRate m_HTF_des = 1000 "Mass flow rate at design point";
      parameter Real gamma = 0.28 "Part of the mass flow going to the recompression directly";
      parameter SI.AngularVelocity[4] choiceN = {75000, 30000, 10000, 3600} * 0.10471975512;
      parameter SI.AngularVelocity N_shaft = choiceN[integer(Modelica.Math.log(P_gro / 10 ^ 6) / Modelica.Math.log(10)) + 2];
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
      parameter Real ratio_m_des = 1 - gamma;
      //Cooler parameters
      parameter SI.ThermodynamicTemperature T_low = 45 + 273.15 "Outlet temperature of the cooler";
      //Exchanger parameters
      parameter SI.ThermodynamicTemperature T_HTF_in_des = 800 + 273.15;
      parameter Integer N_exch = 5;
      //Financial analysis
      parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
      parameter FI.Money C_LTR(fixed = false) "cost of the low temperature heat recuperator";
      parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
      parameter FI.Money C_mainCompressor(fixed = false) "cost of the main compressor";
      parameter FI.Money C_reCompressor(fixed = false) "cost of the re compressor";
      parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
      parameter FI.Money C_generator(fixed = false) "cost of the generator";
      parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
      parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
      parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
      //Results
      SI.Efficiency eta_cycle;
      SI.Energy E_net(final start = 0, fixed = true, displayUnit = "MW.h");
      Boolean m_sup "Disconnect the production of electricity when the outlet pressure of the turbine is close to the critical pressure";
      //Components instanciation
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve HTR(N_q = N_HTR, P_nom_des = P_gro, ratio_m_des = 1) annotation(
        Placement(visible = true, transformation(origin = {12, -22}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.CompressorOnShaft mainCompressor(eta_design = eta_comp_main, N_design = N_shaft, P_nom_des = P_gro, p_high_des = p_high) annotation(
        Placement(visible = true, transformation(origin = {-74, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Cooler cooler(T_low = T_low, P_nom_des = P_gro, T_amb_des = T_amb_des) annotation(
        Placement(visible = true, transformation(origin = {-78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Turbine turbine(PR = PR, N_shaft = N_shaft, eta_design = eta_turb) annotation(
        Placement(visible = true, transformation(origin = {66, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Exchanger exchanger(redeclare package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph, P_nom_des = P_gro, T_out_CO2_des = T_high, N_exch = N_exch, ratio_m_des = 1) annotation(
        Placement(visible = true, transformation(extent = {{34, -8}, {54, 12}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.CompressorOnShaft reCompressor(N_design = N_shaft, P_nom_des = P_gro, p_high_des = p_high) annotation(
        Placement(visible = true, transformation(origin = {-54, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve LTR(N_q = N_LTR, P_nom_des = P_gro, ratio_m_des = 1 - gamma) annotation(
        Placement(visible = true, transformation(origin = {-42, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.FlowMixer mixer annotation(
        Placement(visible = true, transformation(origin = {-20, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.FlowSplitter splitter(gamma = gamma) annotation(
        Placement(visible = true, transformation(origin = {-62, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      parameter MedRec.ThermodynamicState state_HTF_in_des = MedRec.setState_pTX(1.0325 * 10 ^ 5, T_HTF_in_des);
      //     SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = 800+273.15, p_out = 10 ^ 5, m_flow = exchanger.m_HTF_des, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
      //        Placement(visible = true, transformation(origin = {-52, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      //      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sink annotation(
      //        Placement(visible = true, transformation(origin = {-50, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Modelica.Blocks.Interfaces.RealInput parasities_internal;
    protected
    initial equation
      exchanger.h_in_HTF_des = MedRec.specificEnthalpy(state_HTF_in_des);
      exchanger.p_in_HTF_des = state_HTF_in_des.p;
      exchanger.m_HTF_des = m_HTF_des;
      P_nom = (-turbine.W_turb_des) - mainCompressor.W_comp_des - reCompressor.W_comp_des - cooler.P_cool_des;
// enthalpy equalities
//main loop
      exchanger.h_in_CO2_des = HTR.h_out_comp_des;
      turbine.h_in_des = exchanger.h_out_CO2_des;
      HTR.h_in_turb_des = turbine.h_out_des;
      LTR.h_in_turb_des = HTR.h_out_turb_des;
      cooler.h_in_des = LTR.h_out_turb_des;
      mainCompressor.h_in_des = cooler.h_out_des;
      LTR.h_in_comp_des = mainCompressor.h_out_des;
// recompression loop
      reCompressor.h_in_des = LTR.h_out_turb_des;
      HTR.h_in_comp_des = ratio_m_des * LTR.h_out_comp_des + (1 - ratio_m_des) * reCompressor.h_out_des;
//pressure equalities
//main loop
      exchanger.p_in_CO2_des = HTR.p_out_comp_des;
      turbine.p_in_des = exchanger.p_out_CO2_des;
      HTR.p_in_turb_des = turbine.p_out_des;
      LTR.p_in_turb_des = HTR.p_out_turb_des;
      cooler.p_in_des = LTR.p_out_turb_des;
      mainCompressor.p_in_des = cooler.p_out_des;
      LTR.p_in_comp_des = mainCompressor.p_out_des;
//recompression loop
      reCompressor.p_in_des = LTR.p_out_turb_des;
      HTR.p_in_comp_des = ratio_m_des * LTR.p_out_comp_des + (1 - ratio_m_des) * reCompressor.p_out_des;
//mass flow equalities
//main loop
//exchanger.m_CO2_des = HTR.m_comp_des;
      turbine.m_des = exchanger.m_CO2_des;
      HTR.m_turb_des = turbine.m_des;
      LTR.m_turb_des = HTR.m_turb_des;
      cooler.m_des = LTR.m_turb_des * ratio_m_des;
      mainCompressor.m_des = cooler.m_des;
      LTR.m_comp_des = mainCompressor.m_des;
//recompression loop
      HTR.m_comp_des = reCompressor.m_des + LTR.m_comp_des;
      reCompressor.m_des = gamma * LTR.m_turb_des;
// Financial Analysis
      C_HTR = if HTR.T_turb_des[N_HTR] >= 550 + 273.15 then 49.45 * HTR.UA_HTR ^ 0.7544 * (1 + 0.02141 * (HTR.T_turb_des[N_HTR] - 550 - 273.15)) else 49.45 * HTR.UA_HTR ^ 0.7544;
      C_LTR = 49.45 * LTR.UA_HTR ^ 0.7544;
      C_turbine = if exchanger.T_CO2_des[2] >= 550 + 273.15 then 406200 * (-turbine.W_turb_des / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (exchanger.T_CO2_des[2] - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb_des / 10 ^ 6) ^ 0.8;
      C_mainCompressor = 1230000 * (mainCompressor.W_comp_des / 10 ^ 6) ^ 0.3392;
      C_reCompressor = 1230000 * (reCompressor.W_comp_des / 10 ^ 6) ^ 0.3392;
      C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
      C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
      C_exchanger = pri_exchanger * exchanger.Q_HX_des * m_HTF_des / 1000;
      C_PB = (C_HTR + C_LTR + C_turbine + C_mainCompressor + C_reCompressor + C_generator + C_cooler + C_exchanger) * 1.05;
// 1.05 corresponds to inflation from 2017, as correlations are in 2017' dollars.
    equation
      connect(m_sup, exchanger.m_sup) annotation(
        Line);
      connect(exchanger.CO2_port_b, turbine.port_a) annotation(
        Line(points = {{50, -4}, {64.2, -4}, {64.2, -23.8}, {62.2, -23.8}, {62.2, -25.8}}, color = {0, 127, 255}));
      connect(exchanger.CO2_port_a, HTR.from_comp_port_b) annotation(
        Line(points = {{37, -4}, {22.8, -4}, {22.8, -13.6}}, color = {0, 127, 255}));
//  connect(fluid_b, exchanger.HTF_port_b) annotation(
//    Line(points = {{-64, -50}, {12, -50}, {12, 10}, {38, 10}, {38, 10}}, color = {0, 127, 255}));
//  connect(fluid_a, exchanger.HTF_port_a) annotation(
//    Line(points = {{-44, 32}, {54, 32}, {54, 10}, {52, 10}, {52, 10}}, color = {0, 127, 255}));
//connect(src.port_b, fluid_a) annotation(
//  Line(points = {{-60, 64}, {-72, 64}, {-72, 30}, {-42, 30}, {-42, 32}, {-44, 32}}, color = {0, 127, 255}));
//connect(sink.port_a, fluid_b) annotation(
//  Line(points = {{-42, -78}, {-64, -78}, {-64, -54}, {-64, -54}, {-64, -50}}, color = {0, 127, 255}));
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
      connect(turbine.port_b, HTR.from_turb_port_a) annotation(
        Line(points = {{72, -36.6}, {48, -36.6}, {48, -34.6}, {22, -34.6}, {22, -34.1}, {22, -34.1}, {22, -31.6}}, color = {0, 127, 255}));
      connect(cooler.port_b, mainCompressor.port_a) annotation(
        Line(points = {{-78, -46}, {-88, -46}, {-88, -6}, {-82, -6}, {-82, -6}}, color = {0, 127, 255}));
      connect(cooler.T_amb, T_amb);
      connect(m_sup, cooler.m_sup);
//  exchanger.HTF_port_b.m_flow=fluid_b.m_flow;
//  -fluid_a.m_flow+exchanger.HTF_port_a.m_flow=0;
//  exchanger.HTF_port_b.p=fluid_b.p;
//  exchanger.HTF_port_a.p=fluid_a.p;
//  exchanger.HTF_port_a.h_outflow=fluid_a.h_outflow;
//  fluid_b.h_outflow=exchanger.HTF_port_b.h_outflow;
//m_sup = true;
      m_sup = exchanger.HTF_port_a.m_flow >= exchanger.m_HTF_des * nu_min;
      exchanger.CO2_port_a.m_flow = exchanger.m_CO2_des;
      eta_cycle = W_net / exchanger.Q_HX;
      der(E_net) = W_net;
      W_net = if m_sup then (-turbine.W_turb) - mainCompressor.W_comp - reCompressor.W_comp - cooler.P_cooling else 0;
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end recompPB;

    model interpolWithCompo "This model calculates a few off-design points as initial equation and after interpolate over them"
      import SI = Modelica.SIunits;
      import CV = Modelica.SIunits.Conversions;
      import FI = SolarTherm.Models.Analysis.Finances;
      extends SolarTherm.Media.CO2.PropCO2;
      extends Icons.PowerBlock;
      input SI.ThermodynamicTemperature T_amb;
      replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
      replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
      //Objective: having a few points in order to interpolate the HTF outlet temperature and the TIT. Should be feasible with for loops in the initial equation.
      parameter SI.Power P_nom = 100 * 10 ^ 6;
      parameter SI.AbsolutePressure p_high = CV.from_bar(250);
      parameter SI.AbsolutePressure p_low = CV.from_bar(80);
      parameter Real PR = p_high / p_low;
      parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
      parameter SI.ThermodynamicTemperature T_HTF_in_des = CV.from_degC(800) "inlet temperature of the HTF";
      parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
      parameter SI.ThermodynamicTemperature T_amb_des = CV.from_degC(20) "ambient temperature at design";
      parameter SI.MassFlowRate m_des(fixed = false, start = P_nom / 10 ^ 5);
      parameter SI.TemperatureDifference pinchHTR = 5;
      parameter SI.TemperatureDifference pinchExch = T_HTF_in_des - T_high;
      parameter SI.Efficiency eta_comp = 0.89;
      parameter SI.Efficiency eta_turb = 0.93;
      parameter Integer N_exch = 5;
      parameter Integer N_HTR = 15;
      parameter Integer N_cooler = 15;
      //Variables for real time calculation
      //Exchanger variables
      Boolean m_sup "indicates if mass flow from tank is superior enough to calculate";
      SI.MassFlowRate m_HTF_bis "mass flow rate used for switching off the PB";
      SI.ThermodynamicTemperature T_HTF_in "Inlet temperature of HTF";
      SI.ThermodynamicTemperature T_HTF_out "Outlet temperature of HTF";
      SI.ThermodynamicTemperature TIT_CO2 "Turbine inlet Temperature of CO2";
      SI.ThermodynamicTemperature T_CO2_in_HX "CO2 temperature at the inlet of the exchanger";
      MedPB.ThermodynamicState[N_exch] state_CO2_HX(each h.start = 10 ^ 6);
      MedRec.ThermodynamicState[N_exch] state_HTF_HX(each h.start = 6 * 10 ^ 5);
      SI.TemperatureDifference[N_exch] deltaT_HX(each start = 20);
      //PowerBlock variables
      Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
        Placement(transformation(extent = {{-54, 22}, {-34, 42}}), iconTransformation(extent = {{-48, 30}, {-40, 38}})));
      Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
        Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
      Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
        Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
      parameter FI.Money C_PB(fixed = false) "Overall cost of the power block";
      //Exchanger parameters. Calculated only at on-design
      parameter SolarTherm.Types.Conductance[N_exch - 1] UA_exch_dis(each fixed = false, each start = P_nom * 0.4 / N_exch);
      parameter SI.HeatFlowRate Q_dis_exch(fixed = false);
      parameter MedPB.ThermodynamicState[N_exch] state_CO2_exch(each p.fixed = false, each h.fixed = false, each h.start = 10 ^ 6);
      parameter MedRec.ThermodynamicState[N_exch] state_HTF_exch(each p.fixed = false, each h.fixed = false, each h.start = 6 * 10 ^ 5);
      parameter SI.TemperatureDifference[N_exch] deltaT_exch(each start = 20, each fixed = false);
      parameter SI.ThermodynamicTemperature T_HTF_out_des(fixed = false);
      //Off-design parameters for calculation
      parameter Integer N_points = 10;
      parameter SI.ThermodynamicTemperature[N_points] T_high_frac(each fixed = false) "TIT points for off-design";
      parameter SI.ThermodynamicTemperature[N_points] T_CO2_in_HX_table(each fixed = false) "Temperature at the inlet of the exchanger, for a TIT given. Used for interpolation";
      parameter SI.Power[N_points] W_out(each fixed = false) "Outlet power at off-design";
      parameter SI.Efficiency eta_cycle_on(fixed = false, start = 0.5) "efficiency of the cycle";
      parameter SI.Efficiency[N_points] eta_cycle_off(each fixed = false, each start = 0.5) "efficiency of the cycle";
      //Financial analysis
      parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
      parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
      parameter FI.Money C_compressor(fixed = false) "cost of the main compressor";
      parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
      parameter FI.Money C_generator(fixed = false) "cost of the generator";
      parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
      parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.Turbine turbine(eta_design = eta_turb, PR = PR, N_points = N_points) annotation(
        Placement(visible = true, transformation(origin = {38, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Cooler cooler annotation(
        Placement(visible = true, transformation(origin = {-52, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickCO2PB.HeatRecuperatorDTAve HTR(N_points = N_points, N_q = N_HTR, pinchRecuperator = pinchHTR, ratio_m_des = 1, P_nom = P_nom) annotation(
        Placement(visible = true, transformation(origin = {-2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      CompressorOnShaft compr(eta_design = eta_comp, PR = PR, N_points = N_points) annotation(
        Placement(visible = true, transformation(origin = {-38, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
//Part I. On-design=first point
//Inlet of the turbine
      turbine.state_in[1] = MedPB.setState_pTX(p_high, T_high);
      m_des = turbine.m_des;
      turbine.m_des = HTR.m_des;
//High temperature recuperator
      HTR.state_turb[1, N_HTR] = turbine.state_out[1];
      cooler.state_cooler[1] = HTR.state_turb[1, 1];
      cooler.m_des = HTR.m_des;
//Cooler. UA (and therefore this part) are only for economical correlation
      cooler.state_out = MedPB.setState_pTX(p_low, T_low);
//Main compressor. Inlet supposed always the same for off-design
      compr.state_in[1] = cooler.state_out;
      compr.m_des = cooler.m_des;
      HTR.state_comp[1, 1] = compr.state_out[1];
//Exchanger. Done only on-design. After, only the TIT is varied.
      state_CO2_exch[N_exch].h = turbine.state_in[1].h;
      state_CO2_exch[1].h = HTR.state_comp[1, N_HTR].h;
      state_HTF_exch[N_exch] = MedRec.setState_pTX(10 ^ 5, T_HTF_in_des);
      for k in 1:N_exch - 1 loop
        Q_dis_exch = state_CO2_exch[k + 1].h - state_CO2_exch[k].h;
        Q_dis_exch = state_HTF_exch[k + 1].h - state_HTF_exch[k].h;
        m_des * Q_dis_exch = UA_exch_dis[k] * (deltaT_exch[k] + deltaT_exch[k + 1]) / 2;
        deltaT_exch[k] = MedRec.temperature(state_HTF_exch[k]) - MedPB.temperature(state_CO2_exch[k]);
        state_CO2_exch[k].p = p_high;
        state_HTF_exch[k].p = 10 ^ 5;
      end for;
      deltaT_exch[N_exch] = MedRec.temperature(state_HTF_exch[N_exch]) - MedPB.temperature(state_CO2_exch[N_exch]);
      T_HTF_out_des = MedRec.temperature(state_HTF_exch[1]);
      state_CO2_exch[N_exch].p = p_high;
      P_nom = (-turbine.W_turb[1]) - compr.W_comp[1];
      eta_cycle_on = P_nom / (m_des * Q_dis_exch * (N_exch - 1));
//Part I.2. Financial analysis
      C_HTR = if MedPB.temperature(turbine.state_out[1]) >= 550 + 273.15 then 49.45 * sum(HTR.UA_dis) ^ 0.7544 * (1 + 0.02141 * (MedPB.temperature(turbine.state_out[1]) - 550 - 273.15)) else 49.45 * sum(HTR.UA_dis) ^ 0.7544;
      C_turbine = if T_high >= 550 + 273.15 then 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (T_high - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8;
      C_compressor = 1230000 * (compr.W_comp[1] / 10 ^ 6) ^ 0.3392;
      C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
      C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
      C_exchanger = pri_exchanger * m_des * Q_dis_exch * (N_exch - 1) / 1000;
      C_PB = (C_HTR + C_turbine + C_compressor + C_generator + C_cooler + C_exchanger) * 1.05;
//Part II. Off-design points. Mass flow rate is kept constant
      for i in 2:N_points + 1 loop
//    //Inlet of the turbine
        T_high_frac[i - 1] = (0.8 + (1.2 - 0.8) * (i - 1) / (N_points - 1)) * (T_high - 273.15) + 273.15;
        turbine.state_in[i] = MedPB.setState_pTX(p_high, T_high_frac[i - 1]);
//Heat recuperator
        HTR.state_turb[i, N_HTR] = turbine.state_out[i];
        compr.state_in[i] = cooler.state_out;
        HTR.state_comp[i, 1] = compr.state_out[i];
        W_out[i - 1] = (-turbine.W_turb[i]) - compr.W_comp[i];
        eta_cycle_off[i - 1] = W_out[i - 1] / (m_des * (turbine.state_in[i].h - HTR.state_comp[i, N_HTR].h));
        T_CO2_in_HX_table[i - 1] = MedPB.temperature(HTR.state_comp[i, N_HTR]);
      end for;
    equation
      connect(cooler.m_sup, m_sup);
      connect(T_amb, cooler.T_amb);
//Exchanger code
      m_sup = fluid_a.m_flow >= 0.5 * m_des;
      m_HTF_bis = if m_sup then fluid_a.m_flow else m_des;
      for k in 1:N_exch - 1 loop
        m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = m_HTF_bis * (state_HTF_HX[k + 1].h - state_HTF_HX[k].h);
        m_des * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = UA_exch_dis[k] * (deltaT_HX[k] + deltaT_HX[k + 1]) / 2;
        deltaT_HX[k] = MedRec.temperature(state_HTF_HX[k]) - MedPB.temperature(state_CO2_HX[k]);
        state_HTF_HX[k].p = state_HTF_HX[k + 1].p;
        state_CO2_HX[k].p = state_CO2_HX[k + 1].p;
      end for;
      deltaT_HX[N_exch] = MedRec.temperature(state_HTF_HX[N_exch]) - MedPB.temperature(state_CO2_HX[N_exch]);
      state_HTF_HX[N_exch] = if m_sup then MedRec.setState_phX(fluid_a.p, inStream(fluid_a.h_outflow)) else state_HTF_exch[N_exch];
//The TIT and T_CO2 at the inlet of the HX are interpolated from previous points
      T_HTF_in = MedRec.temperature(state_HTF_HX[N_exch]);
      T_HTF_out = MedRec.temperature(state_HTF_HX[1]);
      TIT_CO2 = MedPB.temperature(state_CO2_HX[N_exch]);
      T_CO2_in_HX = Modelica.Math.Vectors.interpolate(T_high_frac, T_CO2_in_HX_table, TIT_CO2);
      state_CO2_HX[1] = MedPB.setState_pTX(p_high, T_CO2_in_HX);
// Power Block calculations
      W_net = if m_sup then max(0, Modelica.Math.Vectors.interpolate(T_high_frac, W_out, TIT_CO2) - cooler.P_cooling) else 0;
//Connectors obligations
      fluid_b.p = fluid_a.p;
      fluid_a.m_flow + fluid_b.m_flow = 0;
      fluid_b.h_outflow = state_HTF_HX[1].h;
      fluid_a.h_outflow = 0;
//shouldn't flow back
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
        __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end interpolWithCompo;
  end Simpler;
  annotation(
      Documentation(info = "<html>
  		<p>This section proposes on-design and off-design calculation of sCO2 cycles. Several off-design points are calculated by varying the turbine's inlet temperature. The exchanger is coded, and allows calculation of the TIT; Modelica is able to solve the entire system. </p>
  		<p>Everything is done in the initial equation section. </p>
  		<p> Simpler proposes a simplified calculation of the different points, by supposing a constant pressure in both high and low pressure of the PB. Complicated proposes the same calculations as DirectDesign. Nevertheless, convergence couldn't be reached when imposing p_out,turb=p_out,comp/PR. Another coding or physical approach is necessary. </p>
  		<p> It was developped from a EES code furnished by the SNL, USA to the CEA, France and from the thesis of J.J. Dyreby, MIT. </p>
  		</html>"));
end QuickCO2PB;