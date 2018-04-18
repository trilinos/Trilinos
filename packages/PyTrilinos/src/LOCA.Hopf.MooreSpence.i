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

%define %loca_hopf_moorespence_docstring
"
PyTrilinos.LOCA.Hopf.MooreSpence is the python interface to namespace
Hopf::MooreSpence of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Hopf.MooreSpence is to provide ***.  The python
version of LOCA.Hopf.MooreSpence supports the following classes:

    * AbstractGroup          - Interface to underlying groups for Hopf point
                               calculations using the Moore-Spence formulation
    * ExtendedGroup          - A group representing the Moore-Spence Hopf
                               equations
    * ExtendedMultiVector    - Multi-vector class to hold solution vectors,
                               Newton vectors, etc.for the Moore-Spence Hopf
                               eqautions
    * ExtendedVector         - Vector class to hold solution vectors, Newton
                               vectors, etc. for Moore-Spence Hopf equations
    * FiniteDifferenceGroup  - Concrete class that provides concrete
                               implementations of the derivative computation
                               methods of the LOCA.Hopf.MooreSpence.-
                               AbstractGroup using first-order finite
                               differencing
    * SolverFactory          - Factory for creating solver objects for solving
                               Moore-Spence Hopf equations
    * SolverStrategy         - Abstract strategy for solving the Moore-Spence
                               Hopf equations
    * SalingerBordering      - Moore-Spence Hopf solver strategy based on
                               'Salinger' bordering.  This is the classic 5-
                               solve Hopf bordering method
"
%enddef

%module(package="PyTrilinos.LOCA.Hopf",
        directors = "1",
        docstring = %loca_hopf_moorespence_docstring) MooreSpence

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

// Trilinos module imports
%import "Teuchos.i"

// Learn about LOCA::Abstract::Iterator::StepStatus enumeration
%import "LOCA_Abstract_Iterator.H"

// Teuchos::RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::MultiContinuation::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::TimeDependent::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::AbstractGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::ExtendedGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::FiniteDifferenceGroup)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SolverStrategy)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SolverFactory)
%teuchos_rcp(LOCA::Hopf::MooreSpence::SalingerBordering)

// Base class support
%pythoncode
%{
import sys, os.path as op
parentDir = op.normpath(op.join(op.dirname(op.abspath(__file__)),".."))
if not parentDir in sys.path: sys.path.append(parentDir)
del sys, op
%}
%import "NOX.Abstract.i"
%import(module="Extended") "LOCA_Extended_MultiAbstractGroup.H"
%import(module="Extended") "LOCA_Extended_MultiVector.H"
%import(module="Extended") "LOCA_Extended_Vector.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="TimeDependent") "LOCA_TimeDependent_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="TurningPoint.MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::Hopf::MooreSpence AbstractGroup class
%include "LOCA_Hopf_MooreSpence_AbstractGroup.H"

// LOCA::Hopf::MooreSpence ExtendedGroup class
%include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"

// LOCA::Hopf::MooreSpence ExtendedMultiVector class
%include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H"

// LOCA::Hopf::MooreSpence ExtendedVector class
%ignore LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency;
%ignore LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam;
%include "LOCA_Hopf_MooreSpence_ExtendedVector.H"

// LOCA::Hopf::MooreSpence FiniteDifferenceGroup class
%include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::Hopf::MooreSpence SolverFactory class
%include "LOCA_Hopf_MooreSpence_SolverFactory.H"

// LOCA::Hopf::MooreSpence SolverStrategy class
%include "LOCA_Hopf_MooreSpence_SolverStrategy.H"

// LOCA::Hopf::MooreSpence SalingerBordering class
%include "LOCA_Hopf_MooreSpence_SalingerBordering.H"
