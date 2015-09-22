// -*- c++ -*-

%module(package="LOCA") Homotopy

%{
// LOCA includes
#include "LOCA_Homotopy_AbstractGroup.H"

// Extra includes due to importing Continuation.i below
#include "LOCA_Continuation_FiniteDifferenceGroup.H"
%}

// Ignore/renames
%ignore operator=;

// Import base class declarations
%import "LOCA_Continuation.i"

// LOCA interface includes
%include "LOCA_Homotopy_AbstractGroup.H"

