//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

// This file should include ALL objects that the user
// is required to instantiate to call LOCA!

// Primary LOCA Objects
#include "NOX.H"
#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Utils.H"

// Bifurcation groups
#include "LOCA_Bifurcation_TPBord_ExtendedVector.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_PitchforkBord_ExtendedVector.H"
#include "LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H"
#include "LOCA_TimeDependent_AbstractGroup.H"
#include "LOCA_Bifurcation_HopfBord_ExtendedVector.H"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
#include "LOCA_Bifurcation_HopfBord_ExtendedGroup.H"

// Homotopy group
#include "LOCA_Homotopy_Group.H"

#include "LOCA_Abstract_Group.H"

// Status tests
#include "LOCA_StatusTest_Wrapper.H"
#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
#include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H"
#include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H"
