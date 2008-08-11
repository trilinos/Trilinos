// -*- c++ -*-

%module(package="LOCA") Bifurcation

%{
// LOCA includes
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.H"
#include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
#include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"
#include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H"
#include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H"
#include "LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H"
#include "LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H"
#include "LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H"

// Extra includes due to importing Continuation.i below
#include "LOCA_Continuation_FiniteDifferenceGroup.H"
#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H"
#include "LOCA_Continuation_StatusTest_ParameterUpdateNorm.H"

#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_FiniteValue.H"
%}

// Ignore/renames
%ignore operator=;
%rename(Print) *::print(ostream& stream, int indent = 0) const;


// Flatten out nested namespaces
%rename(TPBordAbstractGroup) LOCA::Bifurcation::TPBord::AbstractGroup;
%rename(TPBordFiniteDifferenceGroup) LOCA::Bifurcation::TPBord::FiniteDifferenceGroup;
%rename(TPBordSingularSolveGroup) LOCA::Bifurcation::TPBord::SingularSolveGroup;
%rename(HopfBordAbstractGroup) LOCA::Bifurcation::HopfBord::AbstractGroup;
%rename(HopfBordFiniteDifferenceGroup) LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup;
%rename(TPBordNullVectorNormWRMS) LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS;
%rename(TPBordParameterUpdateNorm) LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm;
%rename(PitchforkBordNullVectorNormWRMS) LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS;
%rename(PitchforkBordParameterUpdateNorm) LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm;
%rename(PitchforkBordSlackUpdateNorm) LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm;

// Import base class declarations
%import "LOCA_Continuation.i"
%import "LOCA_TimeDependent.i"

// LOCA interface includes
%include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
%include "LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.H"
%include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"
%include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
%include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"
%include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H"
%include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H"
%include "LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H"
%include "LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H"
%include "LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H"

