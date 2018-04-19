// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
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

%define %nox_import_code
"
from . import ___init__
"
%enddef

%module(package      = "PyTrilinos.NOX",
	autodoc      = "1",
        moduleimport = %nox_import_code,
        docstring    = %nox_docstring) __init__

%{
// System include files
#include <sstream>

// PyTrilinos configuration
#include "PyTrilinos_config.h"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// NOX include files
#include "PyTrilinos_NOX_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// General ignore directives
%ignore operator<<;
%ignore *::operator=;

// Auto-documentation feature
%feature("autodoc", "1");

// Include NOX documentation
%include "NOX_dox.i"

// SWIG library include files
%include "stl.i"

// Trilinos interface import
%import "Teuchos.i"

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

# Import sys module
import sys

# Abstract, Solver, and StatusTest namespaces
__all__ = ['Abstract', 'Solver', 'StatusTest']
from . import Abstract
from . import Solver
from . import StatusTest
%}

// NOX.Epetra namespace
#ifdef HAVE_NOX_EPETRA
%pythoncode
%{

# Epetra namespace
__all__.append('Epetra')
from . import Epetra
%}
#endif

// NOX.PETSc namespace
#ifdef HAVE_NOX_PETSC
%pythoncode
%{

# PETSc namespace
__all__.append('PETSc')
from . import PETSc
sys.modules["PyTrilinos.NOX.PETSc.___init__"] = sys.modules["___init__"]
del sys.modules["___init__"]
%}
#endif
