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

#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"

namespace Rythmos {

enum EStepStatus {
  STEP_STATUS_UNINITIALIZED ///< Stepper is uninitialized
  ,STEP_STATUS_CONVERGED      ///< Nonlinear solver converged and local error test passed
  ,STEP_STATUS_UNKNOWN       ///< Status is unknown
};

inline
const char* toString(const EStepStatus stepStatus)
{
  switch(stepStatus) {
    case STEP_STATUS_UNINITIALIZED: return "STEP_STATUS_UNINITIALIZED";
    case STEP_STATUS_CONVERGED:     return "STEP_STATUS_CONVERGED";
    case STEP_STATUS_UNKNOWN:       return "STEP_STATUS_UNKNOWN";
    default: TEST_FOR_EXCEPT(true);
  }
  return ""; // Never be called!
}

enum EStepLETStatus {
  STEP_LET_STATUS_PASSED     ///< The local truncation error test passed
  ,STEP_LET_STATUS_FAILED    ///< The local truncation error test failed
  ,STEP_LET_STATUS_UNKNOWN   ///< Any local truncation error test was not evaluated
};


inline
const char* toString(const EStepLETStatus stepLETStatus)
{
  switch(stepLETStatus) {
    case STEP_LET_STATUS_PASSED:  return "STEP_LET_STATUS_PASSED";
    case STEP_LET_STATUS_FAILED:  return "STEP_LET_STATUS_FAILED";
    case STEP_LET_STATUS_UNKNOWN: return "STEP_LET_STATUS_UNKNOWN";
    default: TEST_FOR_EXCEPT(true);
  }
  return ""; // Never be called!
}

template<class Scalar>
struct StepStatus {
  std::string message;
  EStepStatus stepStatus;
  EStepLETStatus stepLETStatus;
  Scalar stepSize;
  int order;
  Scalar time;
  Scalar stepLETValue; 
  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > solution;
  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > solutionDot;
  Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > residual;
  Teuchos::RefCountPtr<const Teuchos::ParameterList> extraParameters;
  StepStatus()
    :stepStatus(STEP_STATUS_UNKNOWN)
     ,stepLETStatus(STEP_LET_STATUS_UNKNOWN)
     ,stepLETValue(Scalar(-Teuchos::ScalarTraits<Scalar>::one()))
    {}
};

template<class Scalar>
std::ostream& operator<<( std::ostream& out_arg, const StepStatus<Scalar> &stepStatus )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  *out
    << "message: \"" << stepStatus.message << "\"" << std::endl
    << "stepStatus = " << toString(stepStatus.stepStatus) << std::endl
    << "stepLETStatus = " << toString(stepStatus.stepLETStatus) << std::endl
    << "stepSize = " << stepStatus.stepSize << std::endl
    << "order = " << stepStatus.order << std::endl
    << "time = " << stepStatus.time << std::endl
    << "stepLETValue = " << stepStatus.stepLETValue << std::endl;
  if (stepStatus.solution == Teuchos::null) {
    *out << "solution = NULL" << std::endl;
  } else {
    *out << "solution = " << stepStatus.solution->description();
  }
  if (stepStatus.solutionDot == Teuchos::null) {
    *out << "solutionDot = NULL" << std::endl;
  } else {
    *out << "solutionDot = " << stepStatus.solutionDot->description();
  }
  if (stepStatus.residual == Teuchos::null) {
    *out << "residual = NULL" << std::endl;
  } else {
    *out << "residual = " << stepStatus.residual->description();
  }
  *out << "extraParameters: ";
  if(stepStatus.extraParameters.get()) {
    *out << "\n";
    stepStatus.extraParameters->print(Teuchos::OSTab(out).o(),1000,true);
  }
  else {
    *out << "NONE" << std::endl;
  }
  return out_arg;
};

} // namespace Rythmos

#endif // Rythmos_STEPPER_SUPPORT_TYPES_H



