within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex;

model test
import SI = Modelica.SIunits;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  extends SolarTherm.Media.CO2.PropCO2;
 
  SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = 800+273.15, p_out = 10 ^ 5, m_flow = 1000, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
        Placement(visible = true, transformation(origin = {-52, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {-50, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Complex.simpleRecupPB PB annotation(
    Placement(visible = true, transformation(origin = {42, 6}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
equation
  connect(src.port_b, PB.fluid_a) annotation(
    Line(points = {{-60, 64}, {-58, 64}, {-58, 18}, {26, 18}, {26, 18}}, color = {0, 127, 255}));
  connect(sink.port_a, PB.fluid_b) annotation(
    Line(points = {{-42, -78}, {20, -78}, {20, -10}, {20, -10}}, color = {0, 127, 255}));
  PB.T_amb=300;
end test;