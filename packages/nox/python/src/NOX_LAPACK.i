// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
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

%module(package="PyTrilinos.NOX") LAPACK

// DISREGARD?
// This swig interface file includes the actual definitions for the LAPACK
// module.  See NOX_LAPACK.swi for a description of why this is done this way

%{
// NOX includes
#include "NOX.H"
#include "NOX_LAPACK.H"
%}

// Ignore/renames
//%rename(Print) *::print() const;
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Abstract::Group::operator=(const NOX::LAPACK::Group&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Abstract::Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::LAPACK::Vector&);
%ignore NOX::Abstract::Vector::operator=(const vector<double>&);
%ignore operator<<(ostream&, const NOX::LAPACK::Vector&);
%ignore *::print() const;

// NOX::Abstract imports
%import "NOX_Abstract_Group.H"
%import "NOX_Abstract_Vector.H"

// Import base class declarations
//%import "NOX_Abstract.i"

// LOCA interface includes
%include "NOX_LAPACK_Vector.H"
%include "NOX_LAPACK_Interface.H"
%include "NOX_LAPACK_Group.H"
