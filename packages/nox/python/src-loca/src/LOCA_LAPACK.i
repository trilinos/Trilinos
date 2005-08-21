// -*- c++ -*-

%module(package="PyTrilinos.LOCA") LAPACK

%{
// LOCA includes
#include "LOCA.H"
#include "LOCA_LAPACK.H"
%}

// Ignore/renames
%rename(Print) *::print() const;
%ignore *::operator=;
%ignore operator<<(ostream& stream, const NOX::LAPACK::Vector& v);

// Import LOCA interface
%import "LOCA_Abstract.i"

// Import NOX_LAPACK headers
%import "NOX_LAPACK_Interface.H"
%import "NOX_LAPACK_Group.H"

// LOCA interface includes
%include "LOCA_LAPACK_Interface.H"
%include "LOCA_LAPACK_Group.H"
