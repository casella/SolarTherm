within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex;

 model recompPB "This model calculates a few off-design points as initial equation and after interpolate over them"
 import SI = Modelica.SIunits;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  extends SolarTherm.Media.CO2.PropCO2;
  extends Icons.PowerBlock;
  input SI.ThermodynamicTemperature T_amb;
  replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
  replaceable package MedRec = SolarTherm.Media.SolidParticles.CarboHSP_ph;
  //Objective: having a few points in order to interpolate the HTF outlet temperature and the TIT. Should be feasible with for loops in the initial equation.

  parameter SI.AbsolutePressure p_high = CV.from_bar(200);
  parameter SI.AbsolutePressure p_low = CV.from_bar(80);
  parameter Real PR = p_high / p_low;
  parameter SI.Power P_gro = 100*10^6 "first guess of power outlet";
  parameter SI.Power P_nom (fixed=false) "Electrical power at design point";
  parameter SI.MassFlowRate m_HTF_des = 1000 "Mass flow rate at design point";
  parameter SI.ThermodynamicTemperature T_high = CV.from_degC(715) "inlet temperature of the turbine";
  parameter SI.ThermodynamicTemperature T_HTF_in_des = CV.from_degC(800) "inlet temperature of the HTF";
  parameter SI.ThermodynamicTemperature T_low = CV.from_degC(45) "Inlet temperature of the compressor";
  parameter SI.ThermodynamicTemperature T_amb_des = CV.from_degC(20) "ambient temperature at design";
  parameter Real nu_min=0.25 "Minimum turbine operation" ;
  parameter Real gamma = 0.28 "Part of the mass flow going to the recompression directly";
  parameter SI.MassFlowRate m_des(fixed = false, start = P_nom / 10 ^ 5);
  parameter SI.TemperatureDifference pinchHTR = 5;
  parameter SI.TemperatureDifference pinchLTR = 5;
  parameter SI.TemperatureDifference pinchExch = T_HTF_in_des - T_high;
  parameter SI.Efficiency eta_comp = 0.89;
  parameter SI.Efficiency eta_turb = 0.93;
  parameter Integer N_exch = 5;
  parameter Integer N_HTR = 15;
  parameter Integer N_LTR = 15;
  parameter Integer N_cooler = 15;
  //Variables for real time calculation
  //Exchanger variables
//  Boolean m_sup "indicates if mass flow from tank is superior enough to calculate";
//  SI.MassFlowRate m_HTF_bis "mass flow rate used for switching off the PB";
//  SI.ThermodynamicTemperature T_HTF_in "Inlet temperature of HTF";
//  SI.ThermodynamicTemperature T_HTF_out "Outlet temperature of HTF";
//  SI.ThermodynamicTemperature TIT_CO2 "Turbine inlet Temperature of CO2";
//  SI.ThermodynamicTemperature T_CO2_in_HX "CO2 temperature at the inlet of the exchanger";
//  MedPB.ThermodynamicState[N_exch] state_CO2_HX(each h.start = 10 ^ 6);
//  MedRec.ThermodynamicState[N_exch] state_HTF_HX(each h.start = 6 * 10 ^ 5);
//  SI.TemperatureDifference[N_exch] deltaT_HX(each start = 20);
//  SI.MassFlowRate m_CO2(start = m_des);
//  //PowerBlock variables
//  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = MedRec) annotation(
//    Placement(transformation(extent = {{-54, 22}, {-34, 42}}), iconTransformation(extent = {{-48, 30}, {-40, 38}})));
//  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = MedRec) annotation(
//    Placement(transformation(extent = {{-74, -60}, {-54, -40}}), iconTransformation(extent = {{-62, -48}, {-54, -40}})));
//  Modelica.Blocks.Interfaces.RealOutput W_net(quantity = "Power", unit = "W", displayUnit = "W") "Net electric power output" annotation(
//    Placement(visible = true, transformation(extent = {{78, -22}, {98, -2}}, rotation = 0), iconTransformation(extent = {{46, -10}, {56, 0}}, rotation = 0)));
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
  parameter SI.MassFlowRate[N_points] m_flow(each fixed = false, each start = m_des) "mass flow rate in the cycle for a given off-design point";
  parameter SI.Power[N_points] W_out(each fixed = false) "Outlet power at off-design";
  parameter SI.Efficiency eta_cycle_on(fixed = false, start = 0.5) "efficiency of the cycle";
  //parameter SI.Efficiency[N_points] eta_cycle_off(each fixed = false, each start = 0.5) "efficiency of the cycle";
  //Financial analysis
  parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
  parameter FI.Money C_LTR(fixed = false) "cost of the low temperature heat recuperator";
  parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
  parameter FI.Money C_compressor(fixed = false) "cost of the main compressor";
  parameter FI.Money C_reCompressor(fixed = false) "cost of the main compressor";
  parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
  parameter FI.Money C_generator(fixed = false) "cost of the generator";
  parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
  parameter FI.Money pri_exchanger = 150 "price of the primary exchanger in $/(kW_th). Objective for next-gen CSP with particles";
  
  SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.Turbine turbine(eta_design = eta_turb, PR = PR, N_points = N_points, p_high = p_high) annotation(
    Placement(visible = true, transformation(origin = {38, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.Cooler cooler annotation(
    Placement(visible = true, transformation(origin = {-52, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.HeatRecuperatorDTAve HTR(N_points = N_points, N_q = N_HTR, pinchRecuperator = pinchHTR, ratio_m_des = 1, P_nom = P_nom) annotation(
    Placement(visible = true, transformation(origin = {10, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.CompressorOnShaft compr(eta_design = eta_comp, PR = PR, N_points = N_points, p_high = p_high) annotation(
    Placement(visible = true, transformation(origin = {-32, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.HeatRecuperatorDTAve LTR(N_points = N_points, N_q = N_LTR, pinchRecuperator = pinchLTR, ratio_m_des = 1-gamma, P_nom = P_nom) annotation(
    Placement(visible = true, transformation(origin = {-16, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 CompressorOnShaft reCompr(eta_design = eta_comp, PR = PR, N_points = N_points, p_high = p_high) annotation(
    Placement(visible = true, transformation(origin = {-12, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

initial equation
//Part I. On-design=first point
//Inlet of the turbine
  turbine.state_in[1] = MedPB.setState_pTX(p_high, T_high);
  m_des = turbine.m_flow[1];
  turbine.m_flow[1] = HTR.m_flow_turb[1];
//High and low temperature recuperator
  HTR.state_turb[1, N_HTR] = turbine.state_out[1];
  HTR.state_turb[1, 1] = LTR.state_turb[1,N_LTR];
  cooler.state_cooler[1] = LTR.state_turb[1,1];
  LTR.m_flow_turb[1]=HTR.m_flow_turb[1];
  cooler.m_des = (1-gamma)*LTR.m_flow_turb[1];
//Cooler. UA (and therefore this part) are only for economical correlation
  cooler.state_out = MedPB.setState_pTX(p_low, T_low);
//Main compressor. Inlet supposed always the same for off-design
  compr.state_in[1] = cooler.state_out;
  compr.m_flow[1] = cooler.m_des;
  LTR.state_comp[1, 1] = compr.state_out[1];
  LTR.m_flow_comp[1] = compr.m_flow[1];
  //re compressor loop.
  reCompr.state_in[1]=LTR.state_turb[1,1];
  reCompr.m_flow[1]=gamma*LTR.m_flow_turb[1];
  HTR.state_comp[1,1].h=gamma*reCompr.state_out[1].h+(1-gamma)*LTR.state_comp[1,N_LTR].h;
  HTR.state_comp[1,1].p=gamma*reCompr.state_out[1].p+(1-gamma)*LTR.state_comp[1,N_LTR].p;
  HTR.m_flow_comp[1]=reCompr.m_flow[1]+LTR.m_flow_comp[1];
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
  m_des=m_HTF_des;
  P_nom = (-turbine.W_turb[1]) - compr.W_comp[1]-reCompr.W_comp[1];
  eta_cycle_on = P_nom / (m_des * Q_dis_exch * (N_exch - 1));
//Part I.2. Financial analysis
  C_HTR = if MedPB.temperature(turbine.state_out[1]) >= 550 + 273.15 then 49.45 * sum(HTR.UA_dis) ^ 0.7544 * (1 + 0.02141 * (MedPB.temperature(turbine.state_out[1]) - 550 - 273.15)) else 49.45 * sum(HTR.UA_dis) ^ 0.7544;
  C_LTR = 49.45 * sum(LTR.UA_dis) ^ 0.7544;
  C_turbine = if T_high >= 550 + 273.15 then 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8 * (1 + 1.137 * 10 ^ (-5) * (T_high - 550 - 273.15) ^ 2) else 406200 * (-turbine.W_turb[1] / 10 ^ 6) ^ 0.8;
  C_compressor = 1230000 * (compr.W_comp[1] / 10 ^ 6) ^ 0.3392;
  C_reCompressor = 1230000 * (reCompr.W_comp[1] / 10 ^ 6) ^ 0.3392;
  C_cooler = 32.88 * cooler.UA_cooler ^ 0.75;
  C_generator = 108900 * (P_nom / 10 ^ 6) ^ 0.5463;
  C_exchanger = pri_exchanger * m_des * Q_dis_exch * (N_exch - 1) / 1000;
  C_PB = (C_HTR + C_turbine + C_compressor + C_generator + C_cooler + C_exchanger+C_LTR+C_reCompressor) * 1.05;
//Part II. Off-design points. Mass flow rate is kept constant
  for i in 2:N_points + 1 loop
//Inlet of the turbine
    T_high_frac[i - 1] = (0.8 + (1.2 - 0.8) * (i - 2) / (N_points - 1)) * (T_high - 273.15) + 273.15;
    turbine.T_in[i - 1] = T_high_frac[i - 1];
    turbine.p_in[i - 1] = p_high;
    turbine.p_out[i - 1] = turbine.p_in[i - 1] / PR;
    m_flow[i - 1] = turbine.m_flow[i];
    
//principal loop
    HTR.state_turb[i, N_HTR] = turbine.state_out[i];
    LTR.state_turb[i,N_LTR] = HTR.state_turb[i,1];
    compr.state_in[i] = cooler.state_out;
    LTR.state_comp[i, 1] = compr.state_out[i];
    
    HTR.m_flow_turb[i] = turbine.m_flow[i];
    LTR.m_flow_turb[i]= HTR.m_flow_turb[i];
    compr.m_flow[i] = (1-gamma)*LTR.m_flow_turb[i];
    LTR.m_flow_comp[i] = compr.m_flow[i];
    

     // recompression loop
    reCompr.state_in[i]=LTR.state_turb[i,1];
    reCompr.m_flow[i]=gamma*LTR.m_flow_turb[i];
    HTR.state_comp[i,1].h=gamma*reCompr.state_out[i].h+(1-gamma)*LTR.state_comp[i,N_LTR].h;
    HTR.state_comp[i,1].p=p_high;
    //HTR.state_comp[i,1].p=gamma*reCompr.state_out[i].p+(1-gamma)*LTR.state_comp[i,N_LTR].p;
    HTR.m_flow_comp[i]=reCompr.m_flow[i]+LTR.m_flow_comp[i];
    W_out[i - 1] = (-turbine.W_turb[i]) - compr.W_comp[i]-reCompr.W_comp[i];
////eta_cycle_off[i - 1] = W_out[i - 1] / (m_flow[i-1] * (turbine.state_in[i].h - HTR.state_comp[i, N_HTR].h));
    T_CO2_in_HX_table[i - 1] = MedPB.temperature(HTR.state_comp[i, N_HTR]);
  end for;
equation

  connect(T_amb, cooler.T_amb);
////Exchanger code
//  m_sup = fluid_a.m_flow >= nu_min * m_des;
//  m_HTF_bis = if m_sup then fluid_a.m_flow else m_des;
//  m_CO2 = if m_sup then Modelica.Math.Vectors.interpolate(T_high_frac, m_flow, TIT_CO2) else m_des;
//  for k in 1:N_exch - 1 loop
//    m_CO2 * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = m_HTF_bis * (state_HTF_HX[k + 1].h - state_HTF_HX[k].h);
//    m_CO2 * (state_CO2_HX[k + 1].h - state_CO2_HX[k].h) = UA_exch_dis[k] * (deltaT_HX[k] + deltaT_HX[k + 1]) / 2;
//    deltaT_HX[k] = MedRec.temperature(state_HTF_HX[k]) - MedPB.temperature(state_CO2_HX[k]);
//    state_HTF_HX[k].p = state_HTF_HX[k + 1].p;
//    state_CO2_HX[k].p = state_CO2_HX[k + 1].p;
//  end for;
//  deltaT_HX[N_exch] = MedRec.temperature(state_HTF_HX[N_exch]) - MedPB.temperature(state_CO2_HX[N_exch]);
//  state_HTF_HX[N_exch] = if m_sup then MedRec.setState_phX(fluid_a.p, inStream(fluid_a.h_outflow)) else state_HTF_exch[N_exch];
////The TIT and T_CO2 at the inlet of the HX are interpolated from previous points
//  T_HTF_in = MedRec.temperature(state_HTF_HX[N_exch]);
//  T_HTF_out = MedRec.temperature(state_HTF_HX[1]);
//  TIT_CO2 = MedPB.temperature(state_CO2_HX[N_exch]);
//  T_CO2_in_HX = Modelica.Math.Vectors.interpolate(T_high_frac, T_CO2_in_HX_table, TIT_CO2);
//  state_CO2_HX[1] = MedPB.setState_pTX(p_high, T_CO2_in_HX);
//// Power Block calculations
//  W_net = if m_sup then max(0, Modelica.Math.Vectors.interpolate(T_high_frac, W_out, TIT_CO2) - cooler.P_cooling) else 0;
////Connectors obligations
//  fluid_b.p = fluid_a.p;
//  fluid_a.m_flow + fluid_b.m_flow = 0;
//  fluid_b.h_outflow = state_HTF_HX[1].h;
//  fluid_a.h_outflow = 0;
//shouldn't flow back
  annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
    end recompPB;