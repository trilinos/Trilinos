// -*- c++ -*-

%module(package="PyTrilinos.LOCA") Epetra

%{
// LOCA includes
#include "LOCA.H"
#include "LOCA_Epetra.H"
%}

// Ignore/renames
%rename(Print) *::print() const;
%ignore *::operator=;
%ignore operator<<(ostream& stream, const NOX::Epetra::Vector& v);

// Import LOCA interface
%import "LOCA_Abstract.i"

// Import NOX_Epetra headers
%import "NOX_Epetra_Interface.H"
%import "NOX_Epetra_Group.H"

// LOCA interface includes
%include "LOCA_Epetra_Interface.H"
%include "LOCA_Epetra_Group.H"
