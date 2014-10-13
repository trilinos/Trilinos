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

%define %loca_timedependent_docstring
"
PyTrilinos.LOCA.TimeDependent is the python interface to namespace
TimeDependent of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.TimeDependent is to provide an abstract group for
time dependent problems with a mass matrix.  The python version of
LOCA.TimeDependent supports the following classes:

    * AbstractGroup  - Interface to underlying groups for time dependent
                       systems
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_timedependent_docstring) TimeDependent

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

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Teuchos support
%import "Teuchos.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)

// Import base class declarations
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"

// LOCA::TimeDependent AbstractGroup class
%include "LOCA_TimeDependent_AbstractGroup.H"
