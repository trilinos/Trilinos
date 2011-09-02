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

%define %nox_docstring
"
PyTrilinos.NOX is the python interface to the Trilinos nonlinear
solver package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX is to provide robust nonlinear solvers for the
problem of finding x such that F(x)=0.  In C++, NOX supports several
namespaces, some of which are sub-modules in python:

    * Abstract          - Base classes for abstract interface to NOX
    * Epetra            - Epetra implementation
    * Epetra.Interface  - Concrete interface for Epetra
    * Solver            - Solver manager class and supporting utilities
    * StatusTest        - Support for customizable stopping criteria

The top-level NOX module provides the following user-level class:

    * Utils  - Various utilities

For an example of usage of all of NOX, please consult the following
script in the example subdirectory of the PyTrilinos package:

    * exNOX_1Dfem.py
"
%enddef

%module(package   = "PyTrilinos.NOX",
	autodoc   = "1",
	docstring = %nox_docstring) __init__

%{
// System includes
#include <sstream>

// PyTrilinos configuration
#include "PyTrilinos_config.h"

// Teuchos include
#include "PyTrilinos_Teuchos_Util.h"

// NOX includes
#include "NOX_Version.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// General ignore directives
%ignore operator<<;
%ignore *::operator=;

// Auto-documentation feature
%feature("autodoc", "1");

// Include NOX documentation
%include "NOX_dox.i"

// SWIG library includes
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"
// Note: Teuchos.i turns off warnings for nested classes, so we do not
// have to do it again.

//////////////////////////////////////
// PyTrilinos configuration support //
//////////////////////////////////////
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%constant bool Have_Epetra = true;
#else
%constant bool Have_Epetra = false;
#endif

/////////////////////////
// NOX Version support //
/////////////////////////
%include "NOX_Version.H"
%pythoncode
%{
__version__ = version().split()[2]
%}

///////////////////////
// NOX Utils support //
///////////////////////
%include "NOX_Utils.i"

// NOX namespace imports
%pythoncode
%{
# Abstract, Solver, and StatusTest namespaces
__all__ = ['Abstract', 'Solver', 'StatusTest']
import Abstract
import Solver
import StatusTest
%}

// NOX.Epetra namespace
#ifdef HAVE_NOX_EPETRA
%pythoncode
%{

# Epetra namespace
__all__.append('Epetra')
%}
#endif
