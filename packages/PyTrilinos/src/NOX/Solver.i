// -*- c++ -*-

%module Solver

%{
// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_Manager.H"

// Namespace flattening
using namespace NOX        ;
using namespace NOX::Solver;
%}

// Ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream &, int) const;

// Rename directives
%rename(StatusTest_Generic) NOX::StatusTest::Generic;
%rename(StatusTest_None   ) NOX::StatusTest::None;

// NOX interface includes
%include "NOX_StatusTest_Generic.H"
%include "NOX_Solver_Generic.H"
%include "NOX_Solver_Manager.H"
