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

%define %loca_abstract_docstring
"
PyTrilinos.LOCA.Abstract is the python interface to namespace Abstract
of the Trilinos package LOCA:

    http://trilinos.sandia.gov/packages/nox

"
%enddef

%module(package      = "PyTrilinos.LOCA",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %loca_abstract_docstring) Abstract

%{
// Teuchos includes
#include "Teuchos_PythonParameter.h"

// LOCA includes
#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_Iterator.H"
#include "LOCA_Abstract_TransposeSolveGroup.H"

// Local includes
#include "NumPyImporter.h"
%}

// Standard exception handling
%include "exception.i"

// Include NOX documentation
// %include "LOCA_dox.i"   // TODO: this file will need to be generated

// General ignore directives
%ignore *::operator=;
%ignore *::operator[];

// Trilinos module imports
%import "Teuchos.i"

// Import LOCA interfaces
//%import "LOCA.Continuation.i"
%import "LOCA.MultiContinuation.i"
//%import "LOCA.Homotopy.i"
//%import "LOCA.TimeDependent.i"
//%import "LOCA.Bifurcation.i"

%rename(TimeDependent_AbstractGroup) LOCA::TimeDependent::AbstractGroup;
%rename(Homotopy_AbstractGroup) LOCA::Homotopy::AbstractGroup;
%include "LOCA_TimeDependent_AbstractGroup.H"
%include "LOCA_Homotopy_AbstractGroup.H"

%import "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
%import "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
%import "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
%import "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"
%import "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%import "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
%import "LOCA_Hopf_MooreSpence_AbstractGroup.H"
%import "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
%import "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"
%import "LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H"

// Teuchos::RCPs typemaps
%teuchos_rcp_typemaps(LOCA::GlobalData)
%teuchos_rcp_typemaps(LOCA::DerivUtils)

%include "LOCA_Abstract_Group.H"
%include "LOCA_Abstract_Iterator.H"
%include "LOCA_Abstract_TransposeSolveGroup.H"
