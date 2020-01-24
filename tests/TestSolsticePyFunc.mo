block TestSolsticePyFunc

	function RunPyFunction
		input Integer argc;
		input String argsv[:];
		output Integer result;
		external result = TestExternalPy_func(argc, argsv)
		annotation(Library="python2.7",
			IncludeDirectory="modelica://SolarTherm/Resources/Include",
			Include="#include \"run_py_func.c\""
			);
	end RunPyFunction;

    parameter Integer result(fixed=false);
	parameter String ppath = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Resources/Include/SolsticePy") "Absolute path to the Python script";
	parameter String pname = "run_solstice" "Python script name";
	parameter String pfunc = "run_simul" "Python functiuon name"; // i.e. c = a* b
	parameter Integer argc = 4 "Number of arguments to be passed to the C function";


initial algorithm
	//y = RunPyFunction(argc,{ppath, pname, pfunc, String(time)});
    result := RunPyFunction(argc,{ppath, pname, pfunc, String(time)});
end TestSolsticePyFunc;


