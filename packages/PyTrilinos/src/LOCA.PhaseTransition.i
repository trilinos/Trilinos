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

%define %loca_phasetransition_docstring
"
PyTrilinos.LOCA.PhaseTransition is the python interface to namespace
PhaseTransition of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.PhaseTransition is to provide groups and vectors
for phase transition bifurcations.  The python version of
LOCA.PhaseTransition supports the following classes:

    * AbstractGroup  - Interface to underlying groups for phase transition
                       calculations
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        docstring = %loca_phasetransition_docstring) PhaseTransition

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
#include "LOCA_Extended_MultiAbstractGroup.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Ignore/renames
%ignore *::operator=;

// Standar exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// PyTrilinos module imports
%import "Teuchos.i"

// Teuchos::RCP handling
%teuchos_rcp(LOCA::MultiContinuation::AbstractGroup)
%teuchos_rcp(LOCA::PhaseTransition::AbstractGroup)

// Import base class declarations
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"

// LOCA::PhasTransition AbstractGroup class
%include "LOCA_PhaseTransition_AbstractGroup.H"
