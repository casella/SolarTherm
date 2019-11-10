within SolarTherm.Models.PowerBlocks.sCO2Cycle;
package DirectDesign
  

annotation(
      Documentation(info = "<html>
  		<p>This section proposes direct-design calculation of sCO2 cycles. The on-design is performed in the initial equation section, and off-design is performed in the equation section. </p>
  		<p> It is therefore necessary to connect cycles between them in the initial section as well, through thermodynamic states at the inlet and outlet. </p>
  		<p> It was developped by the CEA (France) from the thesis of J.J. Dyreby, MIT. It still needs validation. </p> 
  		<p> Calculation times are a bit heavy; this should be used for validation of the optimization performed with QuickCO2PB </p>
  		</html>"));
end DirectDesign;