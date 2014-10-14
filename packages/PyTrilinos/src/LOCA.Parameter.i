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

%define %loca_parameter_docstring
"
PyTrilinos.LOCA.Parameter is the python interface to namespace
Parameter of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Parameter is to provide a centralized library for
setting/retrieving numerical parameter values in application codes.
The python version of LOCA.Parameter supports the following classes:

    * Library  - Class to provide a centralized library for setting/retrieving
                 numerical parameter values in application codes
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        docstring = %loca_parameter_docstring) Parameter

%{
// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Ignore/renames
%ignore *::operator=;
%ignore *::operator[];

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Teuchos support
%import "Teuchos.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::Parameter::Library)
%teuchos_rcp(LOCA::Parameter::LibraryT)
%teuchos_rcp(LOCA::Parameter::Vector)

// Import base class declarations

// LOCA::Parameter Library class
%include "LOCA_Parameter_Library.H"

// LOCA::Parameter LibraryT implementations
%include "LOCA_Parameter_LibraryT.H"

// LOCA::Parameter Vector class
%include "LOCA_Parameter_Vector.H"
