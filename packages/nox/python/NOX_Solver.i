// -*- c++ -*-

%module(package="NOX") Solver

%{
// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_Manager.H"
%}

// Ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream& stream, int indent = 0) const;

// Rename directives
%rename(StatusTest_Generic) NOX::StatusTest::Generic;
%rename(StatusTest_None   ) NOX::StatusTest::None;

// NOX interface includes
%include "NOX_StatusTest_Generic.H"
%include "NOX_Solver_Generic.H"
%include "NOX_Solver_Manager.H"
