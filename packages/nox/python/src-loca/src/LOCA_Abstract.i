// -*- c++ -*-

%module(package="PyTrilinos.LOCA") Abstract

%{
// LOCA includes
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_Abstract_Group.H"
%}

// Ignore/renames
%ignore *::operator=;

// Import LOCA interfaces
%import "LOCA_Continuation.i"
%import "LOCA_MultiContinuation.i"
%import "LOCA_Homotopy.i"
%import "LOCA_TimeDependent.i"
%import "LOCA_Bifurcation.i"

// Import NOX headers
%import "NOX_Abstract_Group.H"
%import "NOX_StatusTest_Generic.H"

// LOCA interface includes
%include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%include "LOCA_Abstract_Group.H"

