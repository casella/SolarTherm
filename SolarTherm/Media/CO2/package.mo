within SolarTherm.Media;

package CO2 "Medium model of CO2: bilinear interpolation"
  extends Modelica.Icons.VariantsPackage;




  annotation(
    Documentation(info = "<html>
    <p> The data was created using CoolProp and its wrapper for Mathematica. The Mathematica code and the .txt created are located in SolarTherm/Data/CO2.  </p>
    <p> The range of use is between (P,T)=(60 bar,20°C) and (P,T)=(300 bar, 800°C) </p>
    <p> Call of the data from Modelica was complicated. A solution was found, only working with OpenModelica (and not Dymola, which doesn't enjoy the constant in CO2_ph). It was only checked on Linux system.</p>
    <p>The implementation is taking part of the code of CombiTable2D; declaration of the table are done partly in PropCO2(the table used in setState function); they are declared as outer. Others are declared as constant in CO2_ph. </p>  
		<p>First version released by S. Kamerling for <b>CEA-LITEN/DTBH/SSETI/LSST (France)</b></p>
		</html>"));
end CO2;

