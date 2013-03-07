// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %loca_extended_docstring
"
PyTrilinos.LOCA.Extended is the python interface to namespace Extended
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Extended is to provide classes that extend
NOX.Abstract classes to handle an arbitrary number of multi-vectors
and scalars.  The python version of LOCA.Extended supports the
following classes:

    * MultiVector         - Implemenatation of the NOX.Abstract.MultiVector
                            class for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * Vector              - Implemenatation of the NOX.Abstract.Vector class
                            for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * MultiAbstractGroup  - LOCA abstract interface for extended groups,
                            derived from the NOX.Abstract.Group, i.e., an
                            abstract interface for 'super' groups that have an
                            underlying group component
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_extended_docstring) Extended

%include "LOCA.Extended_Content.i"
