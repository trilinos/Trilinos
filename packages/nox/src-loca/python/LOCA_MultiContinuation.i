// -*- c++ -*-

%module(package="LOCA") MultiContinuation

%{
// LOCA includes
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

// Extra includes due to importing Continuation.i below
#include "LOCA_Continuation_FiniteDifferenceGroup.H"
%}

// Ignore/renames
%ignore operator=;

// Import base class declarations
%import "LOCA_Continuation.i"

// LOCA interface includes
%include "LOCA_MultiContinuation_AbstractGroup.H"
%include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

