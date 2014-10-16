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

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_NumPyVector.hpp"
%}
#endif

// Exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

%import "Teuchos.i"

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
