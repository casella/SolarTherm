within SolarTherm.Models.PowerBlocks;
package sCO2Cycle
  extends Icons.PowerPackage;

annotation(
      Documentation(info = "<html>
  		<p>This section proposes models of sCO2 cycles. <u>No pressure drop is implemented</u> </p>
  		<p> <ul>
  		<li>OnDesign: on-design calculation of sCO2 power blocks </li>
  		<li>OffDesign: off-design calculation, with values from OnDesign to integrate. Useful to debug. </li>
  		<li> DirectDesign: on-design is performed in the initial equation section, off-design in the equation. </li>
  		<li> QuickCO2PB: on-design and some off-design points are calculated in initial section. The exchanger is implemented in the equation section, and other results comes from interpolation. </li>
  		</ul>
  		<p> TODO: implement pressure drops </p>
  		<p> TODO: find a way to converge with the closure equation p_out, turb=p_out,comp/PR. Actually implemented: p_in,turb=p_high. </p>
  		<p> It was developped from a EES code furnished by the SNL, USA and ANU, Australia to the CEA, France, and from the thesis of J.J. Dyreby, MIT. </p>
  		</html>"));
end sCO2Cycle;