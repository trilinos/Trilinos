// -*- c++ -*-

%module(package="LOCA") LAPACK

%{
// LOCA includes
#include "LOCA.H"
#include "LOCA_LAPACK.H"
%}

// Ignore/renames
%rename(Print) *::print() const;
%ignore operator=;
%ignore operator<<(ostream& stream, const NOX::LAPACK::Vector& v);

// Import base class declarations
%import "LOCA_Abstract.i"
%import "NOX_LAPACK_import.i"

// LOCA interface includes
%include "LOCA_LAPACK_Interface.H"
%include "LOCA_LAPACK_Group.H"
