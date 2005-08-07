// -*- c++ -*-

%module(package="PyTrilinos.NOX") Abstract

%{
// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
%}

// Ignore directives
%ignore *::print() const;
%ignore *::operator=;

// Auto-documentation feature
%feature("autodoc", "1");

// NOX interface includes
%include "NOX_Abstract_Group.H"
%include "NOX_Abstract_Vector.H"
