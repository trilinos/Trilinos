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

%define %loca_turningpoint_minimallyaugmented_docstring
"
PyTrilinos.LOCA.TurningPoint.MinimallyAugmented is the python
interface to namespace TurningPoint::MinimallyAugmented of the
Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.TurningPoint.MinimallyAugmented is to provide
groups and vectors for locating turning point bifurcations using the
minimally augmented turning point formulation.  The python version of
LOCA.TurningPoint.MinimallyAugmented supports the following classes:

    * AbstractGroup          - Interface to underlying groups for turning point
                               calculations using the minimally augmented
                               formulation
    * FiniteDifferenceGroup  - Concrete class that provides concrete
                               implementations of the derivative computation
                               methods of the LOCA.TurningPoint.Minimally-
                               Augmented.AbstractGroup using first-order finite
                               differencing
"
%enddef

%module(package   = "PyTrilinos.LOCA.TurningPoint",
        directors = "1",
        docstring = %loca_turningpoint_minimallyaugmented_docstring) MinimallyAugmented

%{
// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

// LOCA includes
#include "LOCA.H"

// Local includes
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
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::AbstractGroup)
%teuchos_rcp(LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup)

// Base class support
%import "NOX.Abstract.i"
%import(module="MultiContinuation") "LOCA_MultiContinuation_AbstractGroup.H"
%import(module="MultiContinuation") "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
%import(module="MooreSpence") "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import(module="MooreSpence") "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"

// LOCA::TurningPoint::MinimallyAugmented AbtractGroup class
%feature("director") LOCA::TurningPoint::MinimallyAugmented;
%include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"

// LOCA::TurningPoint::MinimallyAugmented FinitDifferenceGroup class
%include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
