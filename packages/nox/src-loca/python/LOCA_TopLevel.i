// -*- c++ -*-

%module(package="LOCA") TopLevel

%{
// LOCA includes
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"
#include "Utils_enums.H"
%}

// Ignore/renames
%ignore operator=;
%ignore operator[];
%ignore operator<<;
%rename(Print) LOCA::ParameterVector::print(ostream& stream) const;

// include std_string to convert between python strings and C++ strings
%include "std_string.i"
using namespace std;

// LOCA interface includes
%include "LOCA_Abstract_Iterator.H"
%include "LOCA_Stepper.H"
%include "LOCA_Parameter_Vector.H"
%include "Utils_enums.H"


