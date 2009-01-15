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

%module(package="PyTrilinos.LOCA") Continuation

%{
// LOCA includes
//#include "LOCA_Continuation_AbstractGroup.H"
//#include "LOCA_Continuation_FiniteDifferenceGroup.H"
#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"

// Local includes
#include "NumPyImporter.h"
%}

// Ignore/renames
%ignore *::operator=;
%rename(Print) *::print(ostream& stream, int indent = 0) const;

// Import base class declarations
%import "NOX.Abstract.i"
%import "NOX_StatusTest_Generic.H"

// LOCA interface includes
//%include "LOCA_Continuation_AbstractGroup.H"
//%include "LOCA_Continuation_FiniteDifferenceGroup.H"
//%include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
//%include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
