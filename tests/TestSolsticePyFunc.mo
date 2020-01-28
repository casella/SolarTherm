block TestSolsticePyFunc
	import SolarTherm.Models.CSP.CRS.HeliostatsField.FileOE;

	function RunPyFunction
		input String ppath;
		input String pname;
		input String pfunc;
		input Integer argc;
        input String varnames[:];
        input Real vars[:];
		output Integer result;
		external result = TestExternalPy_func(ppath, pname, pfunc, argc, varnames, vars)
		annotation(Library="python2.7",
			IncludeDirectory="modelica://SolarTherm/Resources/Include",
			Include="#include \"run_py_func.c\""
			);
	end RunPyFunction;

    parameter Integer result(fixed=false);
	parameter String ppath = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Resources/Include/SolsticePy") "Absolute path to the Python script";
	parameter String pname = "run_solstice" "Name of the Python script";
	parameter String pfunc = "run_simul" "Name of the Python functiuon"; 

	parameter Integer argc = 6 "Number of variables to be passed to the C function";

    parameter Integer n_heliostat=20 "Number of heliostats";
    parameter Real H_tower=200 "Tower height";
    parameter Real W_heliostat=10 "Height of a heliostat";
    parameter Real H_heliostat=10 "Width of a heliostat";
    parameter Real rec_h=24.789 "height of the receiver";
    parameter Real rec_w=24.789 "width of the receiver";


    FileOE oeff(nelem=1, file="/home/yewang/solartherm-solstice/tests/results/num_rays_100000/test_solstice.motab");

    //parameter String opt_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Optics/g3p3_opt_eff_1_azim_sud.motab");
    //parameter Solar_angles angles = Solar_angles.ele_azi "Angles used in the lookup table file";

initial algorithm
    result := RunPyFunction(ppath, pname, pfunc, argc, {"n_heliostat","H_tower", "W_heliostat", "H_heliostat", "rec_h", "rec_w"}, {n_heliostat, H_tower, W_heliostat, H_heliostat, rec_h, rec_w});

initial equation
	oeff.wbus.alt = 0.0;
	oeff.wbus.azi = 0.0; 

equation
	  when time >= 3 then
        if result==999 then
		    oeff.wbus.alt = 67.5;
		    oeff.wbus.azi = 0.0;
        else
            oeff.wbus.alt=0.0;
            oeff.wbus.azi=0.0;
        end if; 
	  end when;
end TestSolsticePyFunc;


