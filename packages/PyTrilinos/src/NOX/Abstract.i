// -*- c++ -*-

%module Abstract

%{
// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"

// Namespace flattening
using namespace NOX          ;
using namespace NOX::Abstract;
%}

// Ignore directives
%ignore *::print(ostream &, int) const;
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Abstract::Vector::operator=(const Epetra_Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Epetra::Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Abstract::Vector&);
%ignore NOX::Abstract::Vector::print() const;

// NOX interface includes
%include "NOX_Abstract_Group.H"
%include "NOX_Abstract_Vector.H"
