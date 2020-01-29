block TestSolsticePyFunc
    import SolarTherm.Types.Solar_angles;
    import SolarTherm.Models.CSP.CRS.HeliostatsField.Optical.Table;
    import SolarTherm.Models.CSP.CRS.HeliostatsField.HeliostatsField;

    import SolarTherm.Models.Sources.SolarModel.Sun;
    import nSI = Modelica.SIunits.Conversions.NonSIunits;

	function RunPyFunction
		input String ppath;
		input String pname;
		input String pfunc;
        input String psave;
		input Integer argc;
        input String varnames[:];
        input Real vars[:];
		output Integer result;
		external result = TestExternalPy_func(ppath, pname, pfunc, psave, argc, varnames, vars)
		annotation(Library="python2.7",
			IncludeDirectory="modelica://SolarTherm/Resources/Include",
			Include="#include \"run_py_func.c\""
			);
	end RunPyFunction;

 

	parameter String ppath = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Resources/Include/SolsticePy") "Absolute path to the Python script";
	parameter String pname = "run_solstice" "Name of the Python script";
	parameter String pfunc = "run_simul" "Name of the Python functiuon"; 
    parameter String psave = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Resources/Include/SolsticePy/result/demo") "the directory for saving the results";    

	parameter Integer argc = 6 "Number of variables to be passed to the C function";

    parameter Integer n_heliostat=20 "Number of heliostats";
    parameter Real H_tower=200 "Tower height";
    parameter Real W_heliostat=10 "Height of a heliostat";
    parameter Real H_heliostat=10 "Width of a heliostat";
    parameter Real rec_h=24.789 "height of the receiver";
    parameter Real rec_w=24.789 "width of the receiver";


    parameter String opt_file = psave+"/OELT_Solstice.motab";
    parameter Solar_angles angles = Solar_angles.ele_azi "Angles used in the lookup table file";
    parameter nSI.Angle_deg lat=23.795 "Latitude (+ve North)" annotation(Dialog(group="System location and year"));

    Integer solstice_res;
    Real opt_eff;

    Sun sun(lat=lat, dni=1000);
    Table table(angles = angles, file = opt_file, hra=sun.solar.hra, dec=sun.solar.dec);

initial algorithm
solstice_res := RunPyFunction(ppath, pname, pfunc, psave, argc, {"n_heliostat","H_tower", "W_heliostat", "H_heliostat", "rec_h", "rec_w"}, {n_heliostat, H_tower, W_heliostat, H_heliostat, rec_h, rec_w});

equation
opt_eff = table.nu;

end TestSolsticePyFunc;


