// -*- c++ -*-

%module(package="PyTrilinos.NOX") TopLevel

%{
// System includes
#include <sstream>

// NOX top-level includes
#include "NOX_Version.H"
#include "Utils_enums.H"
%}

// Auto-documentation feature
%feature("autodoc", "1");

// SWIG library includes
%include "std_string.i"

// NOX top-level interface includes
using namespace std;
%include "NOX_Version.H"
%include "Utils_enums.H"

// Python code for the NOX module
%pythoncode %{

__version__ = version().split()[2]

%}
