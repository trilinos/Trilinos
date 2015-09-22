// -*- c++ -*-

%module(package="LOCA") Continuation

%{
// LOCA includes
#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_Continuation_FiniteDifferenceGroup.H"
#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"

// Extra includes due to importing Abstract.i, StatusTest.i below
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"
%}

// Ignore/renames
%ignore operator=;
%rename(Print) *::print(ostream& stream, int indent = 0) const;

// Import base class declarations
%import "NOX_Abstract.i"
%import "NOX_StatusTest.i"

// LOCA interface includes
%include "LOCA_Continuation_AbstractGroup.H"
%include "LOCA_Continuation_FiniteDifferenceGroup.H"
%include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
%include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
