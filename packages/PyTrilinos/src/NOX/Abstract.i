// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

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
