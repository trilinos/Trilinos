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

%module(package="PyTrilinos.LOCA") TimeDependent

%{
// Teuchos includes
#include "PyTrilinos_Teuchos_Util.h"

// LOCA includes
#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_BorderedSystem_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"

#include "LOCA_TimeDependent_AbstractGroup.H"

// Extra includes due to importing Continuation.i below
#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"
%}

// Ignore/renames
%ignore *::operator=;

// Import base class declarations
//%import "LOCA.Continuation.i"
%import "LOCA.MultiContinuation.i"

// LOCA interface includes
%include "LOCA_TimeDependent_AbstractGroup.H"

