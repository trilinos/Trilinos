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

%define %loca_turningpoint_moorespence_docstring
"
PyTrilinos.LOCA.TurningPoint.MooreSpence is the python interface to
namespace TurningPoint::MooreSpence of the Trilinos continuation
algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.TurningPoint.MooreSpence is to provide groups and
vectors for locating turning point bifurcations using the Moore-
Spence formulation.  The python version of
LOCA.TurningPoint.MooreSpence supports the following classes:

    * AbstractGroup          - Interface to underlying groups for turning point
                               calculations using the Moore-Spence formulation
    * FiniteDifferenceGroup  - Concrete class that provides concrete
                               implementations of the derivative computation
                               methods of the LOCA.TurningPoint.MooreSpence.-
                               AbstractGroup using first-order finite
                               differencing
    * SolverFactory          - Factory for creating solver objects for solving
                               Moore-Spence turning point equations
"
%enddef

%module(package   = "PyTrilinos.LOCA.TurningPoint",
        directors = "1",
        docstring = %loca_turningpoint_moorespence_docstring) MooreSpence

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

// LOCA include files
#include "PyTrilinos_LOCA_Headers.hpp"

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;
%ignore operator=;

// Trilinos module imports
%import "Teuchos.i"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::SolverFactory)

// Base class support
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import "LOCA.Parameter.i"

// LOCA::TurningPoint::MooreSpence AbstractGroup class
%include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"

// LOCA::TurningPoint::MooreSpence FiniteDifferenceGroup class
%include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::TurningPoint::MooreSpence SolverFactory class
%include "LOCA_TurningPoint_MooreSpence_SolverFactory.H"
