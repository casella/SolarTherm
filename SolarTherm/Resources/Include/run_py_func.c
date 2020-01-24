#ifndef ST_EXT_PY_FUNC
#define ST_EXT_PY_FUNC

#include <python2.7/Python.h>
#include <stdio.h>

int TestExternalPy_func(int argc, const char *argv[]);

int TestExternalPy_func(int argc, const char *argv[])
{
    int output;
    PyObject *pName, *pModule, *pFunc;
    PyObject *pArgs, *pValue;
    int i;

    if (argc < 3) {
        fprintf(stderr,"Usage: call pythonfile funcname [args]\n");
        return 1;
    }


    Py_Initialize(); /*  Initialize Interpreter  */

    // add the path of the Python function file to the system path
    PyObject *sys_path = PySys_GetObject("path");
    PyList_Append(sys_path, PyString_FromString((char *)argv[0]));

    // name of the Python file
    pName = PyString_FromString(argv[1]);
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, argv[2]);
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(argc - 3);
            for (i = 0; i < argc - 3; ++i) {
                pValue = PyFloat_FromDouble(atof(argv[i + 3])); //atof() converts string to float
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return 1;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            output = PyFloat_AsDouble(pValue);
            fprintf(stderr,"*************\n");
            fprintf(stderr,"Solstice\n");
            printf("The result is: %d\n", output);

            Py_DECREF(pArgs);
            if (pValue != NULL) {
                //printf("Result of call: %ld\n", PyInt_AsLong(pValue));
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
        return 1;
    }
    Py_Finalize();

    return output;
}

#endif
