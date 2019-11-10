within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign;
package Complex "DirectDesign model implemented in initial equation, solved with interpolation."
  extends Icons.PowerPackage;

annotation(
      Documentation(info = "<html>
  		<p>This section proposes models of sCO2 cycles. <u>No pressure drop is implemented</u> </p>
		<p> TODO: implement pressure drops </p>
  		<p> TODO: find a way to converge with the closure equation p_out, turb=p_out,comp. Actually implemented: p_in,turb=p_high. </p>
  		<p> It was developped from a EES code furnished by the SNL, USA and ANU, Australia to the CEA, France, and from the thesis of J.J. Dyreby, MIT. </p>
  		</html>"));
end Complex;