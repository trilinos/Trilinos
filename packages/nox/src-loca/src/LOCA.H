//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

// This file should include ALL objects that the user
// is required to instantiate to call LOCA!

#include "Teuchos_RefCountPtr.hpp"
#include "NOX.H"

// Primary LOCA Objects
#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Stepper.H"
#include "LOCA_Parameter_Vector.H"
#include "LOCA_Parameter_Library.H"

#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ConstraintInterface.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

#include "LOCA_Abstract_TransposeSolveGroup.H"

// Bifurcation groups
#include "LOCA_TimeDependent_AbstractGroup.H"

// Homotopy group
#include "LOCA_Homotopy_Group.H"

#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_Factory.H"

// For Anasazi::MsgType
#ifdef HAVE_LOCA_ANASAZI
#include "AnasaziTypes.hpp"
#endif

// Status tests
// #include "LOCA_StatusTest_Wrapper.H"
// #include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
// #include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
// #include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H"
// #include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H"
// #include "LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H"
// #include "LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H"
// #include "LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H"
