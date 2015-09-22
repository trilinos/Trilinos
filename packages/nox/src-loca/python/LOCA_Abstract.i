// -*- c++ -*-

%module(package="LOCA") Abstract

%{
// LOCA includes
#include "LOCA_Abstract_Group.H"
%}

// Ignore/renames
%ignore operator=;

// Import base class declarations
%import "LOCA_Continuation.i"
%import "LOCA_MultiContinuation.i"
%import "LOCA_Homotopy.i"
%import "LOCA_TimeDependent.i"
%import "LOCA_Bifurcation.i"

// LOCA interface includes
%include "LOCA_Abstract_Group.H"

