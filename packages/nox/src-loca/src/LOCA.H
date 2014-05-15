//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
