// -*- c++ -*-

%module(package="PyTrilinos.LOCA") Chan

%{
// LOCA includes
#include "ChanProblemInterface.H"
%}

// Ignore/renames
%rename(ProblemInterface) ChanProblemInterface;

// Import base class declarations
%import "LAPACK.i"

// LOCA interface includes
%include "ChanProblemInterface.H"

