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

%define %loca_borderedsystem_docstring
"
PyTrilinos.LOCA.BorderedSystem is the python interface to namespace
BorderedSystem of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.BorderedSystem is to provide an interface for
groups that are bordered systems.  The python version of
LOCA.BorderedSystem supports the following classes:

    * AbstractGroup  - Interface for groups that are bordered systems
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_borderedsystem_docstring) BorderedSystem

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
%import "NOX.Abstract.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::BorderedSystem::AbstractGroup)

// LOCA::BorderedSystem AbstractgGroup class
//%feature("director") LOCA::BorderedSystem::AbstractGroup;
//#define result loca_result
%include "LOCA_BorderedSystem_AbstractGroup.H"
//#undef result
