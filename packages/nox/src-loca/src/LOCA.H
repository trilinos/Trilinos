// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This file should include ALL objects that the user
// is required to instantiate to call LOCA!

#include "Teuchos_RCP.hpp"
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
#include "LOCA_Homotopy_DeflatedGroup.H"

#include "LOCA_Abstract_Group.H"
#include "LOCA_Abstract_Factory.H"

// For Anasazi::MsgType
#ifdef HAVE_LOCA_ANASAZI
#include "AnasaziTypes.hpp"
#endif

// Status tests
#include "LOCA_StatusTest_Wrapper.H"
// #include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"
// #include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
// #include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H"
// #include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H"
// #include "LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H"
// #include "LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H"
// #include "LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H"
