// -*- c++ -*-

%module StatusTest

%{
// System includes
#include <sstream>

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"

// Namespace flattening
using namespace NOX            ;
using namespace NOX::StatusTest;
%}

// Ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream &, int) const;

// Rename directives
%rename(StatusTest_None) NOX::StatusTest::None;

// SWIG library includes
%include "std_string.i"

// NOX::StatusTest interface includes
%include "NOX_StatusTest_Generic.H"
%include "NOX_StatusTest_Combo.H"
%include "NOX_StatusTest_NormF.H"
%include "NOX_StatusTest_NormUpdate.H"
%include "NOX_StatusTest_NormWRMS.H"
%include "NOX_StatusTest_MaxIters.H"
%include "NOX_StatusTest_Stagnation.H"
%include "NOX_StatusTest_FiniteValue.H"

// Extensions
%extend NOX::StatusTest::Generic {
  using namespace std;
  string __str__() {
    stringstream os;
    self->print(os);                  // Put the output in os
    string s = os.str();              // Extract the string from os
    return s.substr(0,s.length()-1);  // Return the string minus trailing \n
  }
}
