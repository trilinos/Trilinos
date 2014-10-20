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

%define %loca_borderedsolver_docstring
"
PyTrilinos.LOCA.BorderedSolver is the python interface to namespace
BorderedSolver of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.BorderedSolver is to provide strategies for
solving bordered systems of equations.  The python version of
LOCA.BorderedSolver supports the following classes:

    * Factory  - Factory for creating BorderedSolver strategy objects
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        docstring = %loca_borderedsolver_docstring) BorderedSolver

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

// Namespace flattening
using Teuchos::RCP;
%}

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::BorderedSolver::Factory)

// LOCA::BorderedSolver Factory class
%include "LOCA_BorderedSolver_Factory.H"
