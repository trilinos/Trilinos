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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %loca_hopf_moorespence_docstring
"
PyTrilinos.LOCA.Hopf.MooreSpence is the python interface to namespace
Hopf::MooreSpence of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Hopf.MooreSpence is to provide ***.  The python
version of LOCA.Hopf.MooreSpence supports the following classes:

    * AbstractGroup          -
    * ExtendedGroup          - 
    * ExtendedMultiVector    -
    * ExtendedVector         -
    * FiniteDifferenceGroup  -
    * SolverFactory          -
    * SolverStrategy         -
    * SalingerBordering      -

Any other notes about the package as a whole. . . .
"
%enddef

%module(package="PyTrilinos.LOCA.Hopf",
        docstring = %loca_hopf_moorespence_docstring) MooreSpence

%include "LOCA.Hopf.MooreSpence_Content.i"
