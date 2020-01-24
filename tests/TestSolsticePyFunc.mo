block TestSolsticePyFunc

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

	parameter Integer argc = 4 "Number of variables to be passed to the C function";

    parameter Integer n_heliostat=20 "Number of heliostats";
    parameter Real H_tower=60 "Tower height";
    parameter Real W_heliostat=9.9 "Height of a heliostat";
    parameter Real H_heliostat=6.6 "Width of a heliostat";


initial algorithm
    result := RunPyFunction(ppath, pname, pfunc, argc, {"n_heliostat","H_tower", "W_heliostat", "H_heliostat"}, {n_heliostat, H_tower, W_heliostat, H_heliostat});
end TestSolsticePyFunc;


