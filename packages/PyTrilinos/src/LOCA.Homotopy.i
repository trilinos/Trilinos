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

%define %loca_homotopy_docstring
"
PyTrilinos.LOCA.Homotopy is the python interface to namespace Homotopy
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Homotopy is to provide groups that allow for
Homotopy to be applied.  The python version of LOCA.Homotopy supports
the following classes:

    * Group          - LOCA's Homotopy algorithm
    * AbstractGroup  - Interface to underlying groups for homotopy calculations
    * DeflatedGroup  - LOCA's deflated Homotopy algorithm
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_homotopy_docstring) Homotopy

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

// Local include files
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// PETSc4Py support
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_PETSC
%include "petsc4py/petsc4py.i"
#endif

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

%import "Teuchos.i"

// Learn about LOCA::Abstract::Iterator::StepStatus enumeration
%import "LOCA_Abstract_Iterator.H"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)
%teuchos_rcp(LOCA::BorderedSystem::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::Group)
%teuchos_rcp(LOCA::Homotopy::AbstractGroup)
%teuchos_rcp(LOCA::Homotopy::DeflatedGroup)

// The LOCA::Homotopy classes derive from base classes in other
// modules, so we import those headers here.
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="Extended")          "LOCA_Extended_MultiAbstractGroup.H"
%import(module="BorderedSystem")    "LOCA_BorderedSystem_AbstractGroup.H"

// LOCA::Homotopy Group class
%include "LOCA_Homotopy_Group.H"

// LOCA::Homotopy AbstractGroup class
%include "LOCA_Homotopy_AbstractGroup.H"

// LOCA::Homotopy DeflatedGroup class
%include "LOCA_Homotopy_DeflatedGroup.H"
