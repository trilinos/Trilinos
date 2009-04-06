//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_STEPPER_SUPPORT_TYPES_H
#define Rythmos_STEPPER_SUPPORT_TYPES_H

#include "Rythmos_Types.hpp"
#include "Thyra_VectorBase.hpp"


namespace Rythmos {

    
/** \brief Step type. */
enum StepSizeType { STEP_TYPE_FIXED, STEP_TYPE_VARIABLE };


/** \brief Convert StepSizeType to string. */
inline
const char* toString( const StepSizeType stepSizeType )
{
  switch(stepSizeType) {
    case STEP_TYPE_FIXED:
      return "STEP_TYPE_FIXED";
    case STEP_TYPE_VARIABLE:
      return "STEP_TYPE_VARIABLE";
#ifdef RYTHMOS_DEBUG
    default:
      TEST_FOR_EXCEPT("Invalid enum value!");
#endif
  }
  return 0; // Should never get here!
}


/** \brief . */
enum EStepStatus {
  STEP_STATUS_UNINITIALIZED ///< Stepper is uninitialized
  ,STEP_STATUS_CONVERGED ///< Nonlinear solver converged and local error test passed
  ,STEP_STATUS_UNKNOWN ///< Status is unknown
};


/** \brief . */
inline
const char* toString(const EStepStatus stepStatus)
{
  switch(stepStatus) {
    case STEP_STATUS_UNINITIALIZED: return "STEP_STATUS_UNINITIALIZED";
    case STEP_STATUS_CONVERGED:     return "STEP_STATUS_CONVERGED";
    case STEP_STATUS_UNKNOWN:       return "STEP_STATUS_UNKNOWN";
#ifdef RYTHMOS_DEBUG
    default: TEST_FOR_EXCEPT(true);
#endif
  }
  return ""; // Never be called!
}


/** \brief . */
enum EStepLETStatus {
  STEP_LET_STATUS_PASSED     ///< The local truncation error test passed
  ,STEP_LET_STATUS_FAILED    ///< The local truncation error test failed
  ,STEP_LET_STATUS_UNKNOWN   ///< No local truncation error test was evaluated
};


/** \brief . */
inline
const char* toString(const EStepLETStatus stepLETStatus)
{
  switch(stepLETStatus) {
    case STEP_LET_STATUS_PASSED:  return "STEP_LET_STATUS_PASSED";
    case STEP_LET_STATUS_FAILED:  return "STEP_LET_STATUS_FAILED";
    case STEP_LET_STATUS_UNKNOWN: return "STEP_LET_STATUS_UNKNOWN";
#ifdef RYTHMOS_DEBUG
    default: TEST_FOR_EXCEPT(true);
#endif
  }
  return ""; // Never be called!
}


/** \brief  . */
enum EBreakPointType {
  BREAK_POINT_TYPE_HARD,
  BREAK_POINT_TYPE_SOFT
};


/** \brief . */
inline
const char* toString(const EBreakPointType breakPointType)
{
  switch(breakPointType) {
    case BREAK_POINT_TYPE_HARD:  return "BREAK_POINT_TYPE_HARD";
    case BREAK_POINT_TYPE_SOFT:  return "BREAK_POINT_TYPE_SOFT";
#ifdef RYTHMOS_DEBUG
    default: TEST_FOR_EXCEPT(true);
#endif
  }
  return ""; // Never be called!
}


/** \brief . */
template<class Scalar>
struct StepStatus {
  /** \brief . */
  std::string message;
  /** \brief . */
  EStepStatus stepStatus;
  /** \brief . */
  EStepLETStatus stepLETStatus;
  /** \brief . */
  Scalar stepSize;
  /** \brief . */
  int order;
  /** \brief . */
  Scalar time;
  /** \brief . */
  Scalar stepLETValue; 
  // 2007/05/21: rabartl: ToDo: Change above stepLetValue to ScalarMag
  // 2007/05/21: rabartl: ToDo: We must define what the Local Error Test (LET)
  // is (i.e. what values go into it's computation, what norms are used etc.).
  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > solution;
  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > solutionDot;
  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> > residual;
  /** \brief . */
  RCP<const Teuchos::ParameterList> extraParameters;
  /** \brief . */
  StepStatus()
    :stepStatus(STEP_STATUS_UNKNOWN)
     ,stepLETStatus(STEP_LET_STATUS_UNKNOWN)
     ,stepLETValue(Scalar(-Teuchos::ScalarTraits<Scalar>::one()))
    {}
};


/** \brief . */
template<class Scalar>
std::ostream& operator<<( std::ostream& out_arg, const StepStatus<Scalar> &stepStatus )
{
  using std::endl;
  RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  *out
    << "message: \"" << stepStatus.message << "\"" << endl
    << "stepStatus = " << toString(stepStatus.stepStatus) << endl
    << "stepLETStatus = " << toString(stepStatus.stepLETStatus) << endl
    << "stepSize = " << stepStatus.stepSize << endl
    << "order = " << stepStatus.order << endl
    << "time = " << stepStatus.time << endl
    << "stepLETValue = " << stepStatus.stepLETValue << endl;
  if (stepStatus.solution == Teuchos::null) {
    *out << "solution = NULL" << endl;
  }
  else {
    *out << "solution = " << stepStatus.solution->description() << endl;
  }
  if (stepStatus.solutionDot == Teuchos::null) {
    *out << "solutionDot = NULL" << endl;
  }
  else {
    *out << "solutionDot = " << stepStatus.solutionDot->description() << endl;
  }
  if (stepStatus.residual == Teuchos::null) {
    *out << "residual = NULL" << endl;
  }
  else {
    *out << "residual = " << stepStatus.residual->description() << endl;
  }
  *out << "extraParameters: ";
  if(stepStatus.extraParameters.get()) {
    *out << "\n";
    stepStatus.extraParameters->print(Teuchos::OSTab(out).o(),1000,true);
  }
  else {
    *out << "NONE" << endl;
  }
  return out_arg;
}

} // namespace Rythmos

#endif // Rythmos_STEPPER_SUPPORT_TYPES_H



