within SolarTherm.Models.PowerBlocks.sCO2Cycle.QuickDesign;
package Simple "Simpler model than directDesign implemented in initial equation, solved with interpolation."
  extends Icons.PowerPackage;

annotation(
      Documentation(info = "<html>
  		<p>This section proposes models of sCO2 cycles. <u>No pressure drop is implemented</u> </p>
		<p> This is thought to be easier to manipulate than the complex one. Pressure are kept constants in off-design calculations.</p>
		<p> TODO: implement pressure drops </p>
  		
  		<p> It was developped from a EES code furnished by the SNL, USA and ANU, Australia to the CEA, France, and from the thesis of J.J. Dyreby, MIT. </p>
  		</html>"));
end Simple;