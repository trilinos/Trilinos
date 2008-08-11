// -*- c++ -*-

%module(package="LOCA") Chan

%{
// LOCA includes
#include "ChanProblemInterface.H"
%}

// Ignore/renames
%rename(ProblemInterface) ChanProblemInterface;

// Import base class declarations
%import "LOCA_LAPACK.i"

// LOCA interface includes
%include "ChanProblemInterface.H"

