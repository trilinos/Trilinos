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

%define %loca_hopf_docstring
"
PyTrilinos.LOCA.Hopf is the python interface to namespace Hopf of the
Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Hopf is to provide groups and vectors for Hopf
bifurcations.  The python version of LOCA.Hopf supports the following
sub-modules:

    * MooreSpence         - Groups and vectors for locating Hopf bifurcations
                            using the Moore-Spence formulation
    * MinimallyAugmented  - Groups and vectors for locating Hopf bifurcations
                            using the minimally augmented Hopf formulation

and classes:

    * ComplexMultiVector  - Multi-vector class to hold two multi-vectors to
                            represent a complex multi-vector
    * ComplexVector       - Vector class to hold two vectors to represent a
                            complex vector
"
%enddef

%module(package   = "PyTrilinos.LOCA.Hopf",
        docstring = %loca_hopf_docstring) __init__

%{
// PyTrilinos include files
#include "PyTrilinos_config.h"
#include "PyTrilinos_LinearProblem.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"
#endif

// NOX-Epetra include files
#ifdef HAVE_NOX_EPETRA
//#include "Epetra_Vector.h"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#endif

// NOX-PETSc include files
#include "NOX_Abstract_Vector.H"
#ifdef HAVE_NOX_PETSC
#include "NOX_Petsc_Vector.H"
#endif

// LOCA include files
#include "PyTrilinos_LOCA_Headers.hpp"
#include "PyTrilinos_LOCA_Hopf_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// PETSc4Py support
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_PETSC
%include "petsc4py/petsc4py.i"
#endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;
%ignore *::operator[];
%ignore operator[];

// Trilinos module imports
%import "Teuchos.i"

// Base class imports
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
if "delete_ComplexMultiVector" not in dir(___init__):
    del ___init__
    from . import ___init__
%}
%import "NOX.Abstract.i"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"

// Import the sub-modules
%pythoncode
%{
__all__ = ['MooreSpence',
           'MinimallyAugmented'
           ]
import MooreSpence
import MinimallyAugmented
%}

// LOCA::Hopf ComplexMultiVector class
%include "LOCA_Hopf_ComplexMultiVector.H"

// LOCA::Hopf ComplexVector class
%include "LOCA_Hopf_ComplexVector.H"
