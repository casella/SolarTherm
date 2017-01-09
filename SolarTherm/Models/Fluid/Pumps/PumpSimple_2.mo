within SolarTherm.Models.Fluid.Pumps;
model PumpSimple_2
  extends SolarTherm.Interfaces.Models.Pump;
  Modelica.Blocks.Interfaces.RealInput m_flow annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-10,78}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={0,86})));
  parameter Real k_loss(unit="J/kg")=0.55e3;
  SI.Power W_loss;
  
  SI.Power Ex_in "Inlet exergy flow rate"; 
  SI.Power Ex_out "Outlet exergy flow rate";
  SI.Power Ex_des "Exergy destrcution rate";
  SI.Energy Ex_loss(final start=0) "Exergy destrcution";
  
  Modelica.Blocks.Interfaces.RealInput T_amb annotation (Placement(
        iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={0,86})));
        
  Medium.ThermodynamicState state_in=Medium.setState_phX(fluid_a.p,inStream(fluid_a.h_outflow));
  Medium.ThermodynamicState state_out=Medium.setState_phX(fluid_b.p,actualStream(fluid_b.h_outflow));
  SI.SpecificEntropy s_in=Medium.specificEntropy(state_in);
  SI.SpecificEntropy s_out=Medium.specificEntropy(state_out);
  Medium.ThermodynamicState state_dead=Medium.setState_pTX(1e5,T_amb);
  SI.SpecificEntropy s_dead=Medium.specificEntropy(state_dead);
  SI.SpecificEnthalpy h_dead=Medium.specificEnthalpy(state_dead);


equation
  fluid_b.m_flow=-m_flow;
  fluid_a.m_flow=m_flow;
  fluid_a.h_outflow=fluid_b.h_outflow;
  fluid_b.h_outflow=inStream(fluid_a.h_outflow);
  fluid_a.Xi_outflow=fluid_b.Xi_outflow;
  fluid_b.Xi_outflow=inStream(fluid_a.Xi_outflow);
  //fluid_a.p=fluid_b.p;
  W_loss=k_loss*m_flow;
  
  Ex_in=fluid_a.m_flow*((inStream(fluid_a.h_outflow)-h_dead)-T_amb*(s_in-s_dead));
  Ex_out=-fluid_b.m_flow*((actualStream(fluid_b.h_outflow)-h_dead)-T_amb*(s_out-s_dead));
  
  Ex_in-Ex_out-Ex_des=0;
  der(Ex_loss)=Ex_des;
  
  annotation (Documentation(revisions="<html>
<ul>
<li>Alberto de la Calle:<br>Released first version. </li>
</ul>
</html>"));
end PumpSimple_2;
