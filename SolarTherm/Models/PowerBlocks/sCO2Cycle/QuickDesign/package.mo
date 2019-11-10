within SolarTherm.Models.PowerBlocks.sCO2Cycle;
package QuickDesign "offer quick calculation for on and off-design through interpolation of the results"
  extends Icons.PowerPackage;

annotation(
      Documentation(info = "<html>
  		<p>This section proposes models of sCO2 cycles. <u>No pressure drop is implemented</u> </p>
		<p>On-design and some off-design points are calculated in initial section. The exchanger is implemented in the equation section, and other results comes from interpolation.</p>
  		<p> <ul>
  		<li> Simple: for off-design analysis, pressures are kept constant, no elliptic equation.  </li>
		<li> Complicated: elliptic equation is implemented. Nevertheless, the cycle is hard to close regardign convergence.</li> 
  		</ul>
  		<p> TODO: implement pressure drops </p>
  		<p> TODO: find a way to converge with the closure equation p_out, turb=p_out,comp. Actually implemented: p_in,turb=p_high. </p>
  		<p> It was developped from a EES code furnished by the SNL, USA and ANU, Australia to the CEA, France, and from the thesis of J.J. Dyreby, MIT. </p>
  		</html>"));
end QuickDesign;