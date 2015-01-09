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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_FIXED_STEP_CONTROL_STRATEGY_DEF_H
#define Rythmos_FIXED_STEP_CONTROL_STRATEGY_DEF_H

#include "Rythmos_FixedStepControlStrategy_decl.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {


// Constructors

template<class Scalar>
void FixedStepControlStrategy<Scalar>::setStepControlState_(
  StepControlStrategyState newState)
{
  if (stepControlState_ == UNINITIALIZED) {
    TEUCHOS_TEST_FOR_EXCEPTION(newState != BEFORE_FIRST_STEP, std::logic_error,
                               "newState = " << toString(newState) << "\n");
  } else if (stepControlState_ == BEFORE_FIRST_STEP) {
    TEUCHOS_TEST_FOR_EXCEPTION(newState != MID_STEP, std::logic_error,
                               "newState = " << toString(newState) << "\n");
  } else if (stepControlState_ == MID_STEP) {
    TEUCHOS_TEST_FOR_EXCEPTION(newState != AFTER_CORRECTION, std::logic_error,
                               "newState = " << toString(newState) << "\n");
  } else if (stepControlState_ == AFTER_CORRECTION) {
    TEUCHOS_TEST_FOR_EXCEPTION(newState != READY_FOR_NEXT_STEP,std::logic_error,
                               "newState = " << toString(newState) << "\n");
  } else if (stepControlState_ == READY_FOR_NEXT_STEP) {
    TEUCHOS_TEST_FOR_EXCEPTION(newState != MID_STEP, std::logic_error,
                               "newState = " << toString(newState) << "\n");
  }
  stepControlState_ = newState;
}

template<class Scalar>
StepControlStrategyState FixedStepControlStrategy<Scalar>::getCurrentState()
{
  return(stepControlState_);
}

template<class Scalar>
FixedStepControlStrategy<Scalar>::FixedStepControlStrategy()
  : stepControlState_(UNINITIALIZED)
{}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::initialize(
  const StepperBase<Scalar>& stepper)
{
  stepControlState_ = UNINITIALIZED;
  // Any other initialization goes here.
  setStepControlState_(BEFORE_FIRST_STEP);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::setRequestedStepSize(
  const StepperBase<Scalar>& stepper,
  const Scalar& stepSize,
  const StepSizeType& stepSizeType)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !((stepControlState_ == UNINITIALIZED) ||
      (stepControlState_ == BEFORE_FIRST_STEP) ||
      (stepControlState_ == READY_FOR_NEXT_STEP) ||
      (stepControlState_ == MID_STEP)), std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::setRequestedStepSize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    (stepSizeType != STEP_TYPE_FIXED),
    std::logic_error, "Error, step size type != STEP_TYPE_FIXED "
    "for FixedStepControlStrategy!\n");

  if (stepControlState_ == UNINITIALIZED) setStepControlData(stepper);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::nextStepSize(
  const StepperBase<Scalar>& stepper, Scalar* stepSize,
  StepSizeType* stepSizeType, int* order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!((stepControlState_ == BEFORE_FIRST_STEP) ||
                               (stepControlState_ == MID_STEP) ||
                               (stepControlState_ == READY_FOR_NEXT_STEP) ),
                               std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::nextStepSize()\n");

  setStepControlState_(MID_STEP);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::setCorrection(
     const StepperBase<Scalar>& stepper
    ,const RCP<const Thyra::VectorBase<Scalar> >& soln
    ,const RCP<const Thyra::VectorBase<Scalar> >& dx
    ,int solveStatus)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != MID_STEP, std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::setCorrection()\n");
  setStepControlState_(AFTER_CORRECTION);
}

template<class Scalar>
bool FixedStepControlStrategy<Scalar>::acceptStep(
  const StepperBase<Scalar>& stepper, Scalar* value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::completeStep()\n");

  return true;
}

template<class Scalar>
AttemptedStepStatusFlag FixedStepControlStrategy<Scalar>::rejectStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::completeStep()\n");

  setStepControlState_(READY_FOR_NEXT_STEP);

  return (REP_ERR_FAIL);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::completeStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FixedStepControlStrategy<Scalar>::completeStep()\n");

  setStepControlState_(READY_FOR_NEXT_STEP);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::describe(
  Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::as;
  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
       (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )    ) {
    out << this->description() << "::describe" << "\n";
  }
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
FixedStepControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }

  return (validPL);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
FixedStepControlStrategy<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
FixedStepControlStrategy<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template<class Scalar>
void FixedStepControlStrategy<Scalar>::setStepControlData(
  const StepperBase<Scalar>& stepper)
{
  if (stepControlState_ == UNINITIALIZED) initialize(stepper);
}

template<class Scalar>
bool FixedStepControlStrategy<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
FixedStepControlStrategy<Scalar>::cloneStepControlStrategyAlgorithm() const
{

  RCP<FixedStepControlStrategy<Scalar> >
    stepControl = rcp(new FixedStepControlStrategy<Scalar>());

  if (!is_null(parameterList_))
    stepControl->setParameterList(parameterList_);

  return stepControl;
}

template<class Scalar>
int FixedStepControlStrategy<Scalar>::getMaxOrder() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepControlState_ == UNINITIALIZED, std::logic_error,
      "Error, attempting to call getMaxOrder before initialization!\n"
      );
  return(0);
}

//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_FIXED_STEP_CONTROL_STRATEGY_INSTANT(SCALAR) \
  template class FixedStepControlStrategy< SCALAR >;


} // namespace Rythmos
#endif // Rythmos_FIXED_STEP_CONTROL_STRATEGY_DEF_H
