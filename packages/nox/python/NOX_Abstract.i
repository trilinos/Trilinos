// -*- c++ -*-

%module(package="NOX") Abstract

%{
// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
%}

// Ignore directives
%ignore *::print() const;
%ignore operator=;

// NOX interface includes
%include "NOX_Abstract_Group.H"
%include "NOX_Abstract_Vector.H"
