within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple;

model test
import SI = Modelica.SIunits;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  extends SolarTherm.Media.CO2.PropCO2;
 
  SolarTherm.Models.PowerBlocks.sCO2Cycle.SourceFlow src(T_out = 800+273.15, p_out = 10 ^ 5, m_flow = 1000, redeclare package MedPB = SolarTherm.Media.SolidParticles.CarboHSP_ph, use_m_parameter = true) annotation(
        Placement(visible = true, transformation(origin = {-52, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      
      SolarTherm.Models.PowerBlocks.sCO2Cycle.SinkFlow sink annotation(
        Placement(visible = true, transformation(origin = {-50, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign.Simple.SimpleRecupPB PB annotation(
    Placement(visible = true, transformation(origin = {46, -8}, extent = {{-58, -58}, {58, 58}}, rotation = 0)));
equation
  connect(sink.port_a, PB.fluid_b) annotation(
    Line(points = {{-42, -78}, {12, -78}, {12, -34}, {12, -34}}, color = {0, 127, 255}));
  connect(src.port_b, PB.fluid_a) annotation(
    Line(points = {{-60, 64}, {20, 64}, {20, 12}, {20, 12}}, color = {0, 127, 255}));
  connect(PB.fluid_b, sink.port_a) annotation(
    Line(points = {{-10, -20}, {-8, -20}, {-8, -78}, {-42, -78}, {-42, -78}}, color = {0, 127, 255}));
  connect(src.port_b, PB.fluid_a) annotation(
    Line(points = {{-60, 64}, {-78, 64}, {-78, 24}, {-2, 24}, {-2, 26}, {-2, 26}}, color = {0, 127, 255}));
  PB.T_amb=300;
end test;