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

%define %loca_multicontinuation_docstring
"
PyTrilinos.LOCA.MultiContinuation is the python interface to namespace
MultiContinuation of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.MultiContinuation is to provide groups and vectors
for multi-parameter continuation.  The python version of
LOCA.MultiContinuation supports the following classes:

    * AbstractGroup            - LOCA abstract interface for continuation,
                                 derived from the NOX.Abstract.Group.  This
                                 abstract class provides the interface
                                 necessary to perform continuation, i.e.,
                                 compute families of solutions to F(x,p) = 0
    * FiniteDifferenceGroup    - Concrete class that provides a concrete
                                 implementation of the computeDfDp() method of
                                 the LOCA.Continuation.AbstractGroup using
                                 first-order finite differencing
    * ConstraintInterface      - Abstract interface for the constraint portion
                                 of a constrained nonlinear system
    * ConstraintInterfaceMVDX  - Abstract interface for the constraint portion
                                 of a constrained nonlinear system for
                                 constraints that support computing a solution
                                 component derivative as a multi-vector
    * ExtendedMultiVector      - MultiVector class to hold solution vectors,
                                 Newton vectors, etc. for continuation equations
    * ExtendedVector           - Vector class to hold solution vectors, Newton
                                 vectors, etc. for continuation equations
    * Factory                  - Factory for creating continuation strategy
                                 objects
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_multicontinuation_docstring) MultiContinuation

%{
// PyTrilinos include files
#include "PyTrilinos_config.h"
#include "PyTrilinos_LinearProblem.hpp"

// Teuchos include files
#include "PyTrilinos_Teuchos_Headers.hpp"

// Epetra include files
#ifdef HAVE_PYTRILINOS_EPETRA
#include "PyTrilinos_Epetra_Headers.hpp"
#endif

// NOX-Epetra include files
#ifdef HAVE_PYTRILINOS_NOX_EPETRA
//#include "Epetra_Vector.h"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#endif

// NOX-PETSc include files
#include "NOX_Abstract_Vector.H"
#ifdef HAVE_PYTRILINOS_NOX_PETSC
#include "NOX_Petsc_Vector.H"
#endif

// LOCA include files
#include "PyTrilinos_LOCA_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// PETSc4Py support
%include "PyTrilinos_config.h"
#ifdef HAVE_PYTRILINOS_NOX_PETSC
%include "petsc4py/petsc4py.i"
#endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

%import "Teuchos.i"

// Learn about LOCA::Abstract::Iterator::StepStatus enumeration
%import "LOCA_Abstract_Iterator.H"

// Teucho::RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterface)
%teuchos_rcp(LOCA::MultiContinuation::ConstraintInterfaceMVDX)
%teuchos_rcp(LOCA::MultiContinuation::Factory)

// Allow import from this directory
%pythoncode
%{
import sys, os.path as op
thisDir = op.dirname(op.abspath(__file__))
if not thisDir in sys.path: sys.path.append(thisDir)
del sys, op
%}

// Import base class declarations
%import "NOX.Abstract.i"
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"

// LOCA::MultiContinuation AbstractGroup class
%include "LOCA_MultiContinuation_AbstractGroup.H"

// LOCA::MultiContinuation FiniteDifferenceGroup class
%include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

// LOCA::MultiContinuation ConstraintInterface class
%feature("director") LOCA::MultiContinuation::ConstraintInterface;
%include "LOCA_MultiContinuation_ConstraintInterface.H"

// LOCA::MultiContinuation ConstraintInterfaceMVDX class
%warnfilter(473) LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%feature("director") LOCA::MultiContinuation::ConstraintInterfaceMVDX;
%include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

// LOCA::MultiContinuation ExtendedMultiVector class
%include "LOCA_MultiContinuation_ExtendedMultiVector.H"

// LOCA::MultiContinuation ExtendedVector class
%include "LOCA_MultiContinuation_ExtendedVector.H"

// LOCA::MultiContinuation Factory class
%include "LOCA_MultiContinuation_Factory.H"
