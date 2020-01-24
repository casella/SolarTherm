block TestSolsticePyFunc

	function RunPyFunction
		input Integer argc;
		input String argsv[:];
        input String varnames[:];
        input Real vars[:];
		output Integer result;
		external result = TestExternalPy_func(argc, argsv, varnames, vars)
		annotation(Library="python2.7",
			IncludeDirectory="modelica://SolarTherm/Resources/Include",
			Include="#include \"run_py_func.c\""
			);
	end RunPyFunction;

    parameter Integer result(fixed=false);
	parameter String ppath = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Resources/Include/SolsticePy") "Absolute path to the Python script";
	parameter String pname = "run_solstice" "Python script name";
	parameter String pfunc = "run_simul" "Python functiuon name"; // i.e. c = a* b
	parameter Integer argc = 4 "Number of variables to be passed to the C function";
    parameter Integer n_heliostat=20 "Number of heliostats";
    parameter Real H_tower=60 "Tower height";
    parameter Real W_heliostat=9.9 "Height of a heliostat";
    parameter Real H_heliostat=6.6 "Width of a heliostat";


initial algorithm
	//y = RunPyFunction(argc,{ppath, pname, pfunc, String(time)});
    result := RunPyFunction(argc,{ppath, pname, pfunc}, {"n_heliostat","H_tower", "W_heliostat", "H_heliostat"}, {n_heliostat, H_tower, W_heliostat, H_heliostat});
end TestSolsticePyFunc;


