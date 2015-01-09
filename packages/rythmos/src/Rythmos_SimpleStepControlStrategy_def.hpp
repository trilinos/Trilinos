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

#ifndef Rythmos_SIMPLE_STEP_CONTROL_STRATEGY_DEF_H
#define Rythmos_SIMPLE_STEP_CONTROL_STRATEGY_DEF_H

#include "Rythmos_SimpleStepControlStrategy_decl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {

// Static members


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::initialStepSizeName_
= "Initial Step Size";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::initialStepSizeDefault_
= std::numeric_limits<Scalar>::min();


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::minStepSizeName_
= "Min Step Size";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::minStepSizeDefault_
= std::numeric_limits<Scalar>::min();


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::maxStepSizeName_
= "Max Step Size";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::maxStepSizeDefault_
= std::numeric_limits<Scalar>::max();


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::stepSizeDecreaseFactorName_
= "Step Size Decrease Factor";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::stepSizeDecreaseFactorDefault_
= 0.5;


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::stepSizeIncreaseFactorName_
= "Step Size Increase Factor";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::stepSizeIncreaseFactorDefault_
= 1.5;


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::maxStepFailuresName_
= "Maximum Number of Step Failures";

template<class Scalar>
const int
SimpleStepControlStrategy<Scalar>::maxStepFailuresDefault_
= 100;


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::dxRelativeToleranceName_
= "Solution Change Relative Tolerance";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::dxRelativeToleranceDefault_
= 1.0e-06;


template<class Scalar>
const std::string
SimpleStepControlStrategy<Scalar>::dxAbsoluteToleranceName_
= "Solution Change Absolute Tolerance";

template<class Scalar>
const double
SimpleStepControlStrategy<Scalar>::dxAbsoluteToleranceDefault_
= 1.0e-12;


// Constructors

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::setStepControlState_(StepControlStrategyState newState)
{
  if (stepControlState_ == UNINITIALIZED) {
    TEUCHOS_TEST_FOR_EXCEPT(newState != BEFORE_FIRST_STEP);
  } else if (stepControlState_ == BEFORE_FIRST_STEP) {
    TEUCHOS_TEST_FOR_EXCEPT(newState != MID_STEP);
  } else if (stepControlState_ == MID_STEP) {
    TEUCHOS_TEST_FOR_EXCEPT(newState != AFTER_CORRECTION);
  } else if (stepControlState_ == AFTER_CORRECTION) {
    TEUCHOS_TEST_FOR_EXCEPT(newState != READY_FOR_NEXT_STEP);
  } else if (stepControlState_ == READY_FOR_NEXT_STEP) {
    TEUCHOS_TEST_FOR_EXCEPT(newState != MID_STEP);
  }
  stepControlState_ = newState;
}

template<class Scalar>
StepControlStrategyState SimpleStepControlStrategy<Scalar>::getCurrentState()
{
  return(stepControlState_);
}

template<class Scalar>
SimpleStepControlStrategy<Scalar>::SimpleStepControlStrategy()
  : stepControlState_(UNINITIALIZED),
    initialStepSize_(initialStepSizeDefault_),
    stepSizeType_(STEP_TYPE_VARIABLE),
    minStepSize_(minStepSizeDefault_),
    maxStepSize_(maxStepSizeDefault_),
    stepSizeIncreaseFactor_(stepSizeIncreaseFactorDefault_),
    stepSizeDecreaseFactor_(stepSizeDecreaseFactorDefault_),
    numStepFailures_(0),
    maxStepFailures_(maxStepFailuresDefault_),
    maxOrder_(1),
    dxRelativeTolerance_(dxRelativeToleranceDefault_),
    dxAbsoluteTolerance_(dxAbsoluteToleranceDefault_),
    solveStatus_(0)
{}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::initialize(
  const StepperBase<Scalar>& stepper)
{
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool doTrace = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH));
  Teuchos::OSTab ostab(out,1,"initialize");

  if (doTrace) {
    *out << "\nEntering " << this->Teuchos::Describable::description()
         << "::initialize()...\n";
  }

  setStepControlState_(BEFORE_FIRST_STEP);

  if (doTrace) {
    *out << "\nLeaving " << this->Teuchos::Describable::description()
         << "::initialize()...\n";
  }
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::setRequestedStepSize(
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
     << ") for SimpleStepControlStrategy<Scalar>::setRequestedStepSize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ((stepSizeType == STEP_TYPE_FIXED) && (stepSize == ST::zero())),
    std::logic_error, "Error, step size type == STEP_TYPE_FIXED, "
    "but requested step size == 0!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ((stepSizeType == STEP_TYPE_FIXED) && (stepSize < minStepSize_)),
    std::logic_error, "Error, step size type == STEP_TYPE_FIXED, "
    "and (stepSize="<<stepSize<<") < (minStepSize="<<minStepSize_<<")!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ((stepSizeType == STEP_TYPE_FIXED) && (stepSize > maxStepSize_)),
    std::logic_error, "Error, step size type == STEP_TYPE_FIXED, "
    "and (stepSize="<<stepSize<<") > (maxStepSize="<<maxStepSize_<<")!\n");

  if (stepControlState_ == UNINITIALIZED) initialize(stepper);
  requestedStepSize_ = stepSize;
  stepSizeType_ = stepSizeType;
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::nextStepSize(
  const StepperBase<Scalar>& stepper, Scalar* stepSize,
  StepSizeType* stepSizeType, int* order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!((stepControlState_ == BEFORE_FIRST_STEP) ||
                               (stepControlState_ == MID_STEP) ||
                               (stepControlState_ == READY_FOR_NEXT_STEP) ),
                               std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for SimpleStepControlStrategy<Scalar>::nextStepSize()\n");

  if (stepControlState_ == BEFORE_FIRST_STEP) {
    if (initialStepSize_ == initialStepSizeDefault_)
      initialStepSize_ = requestedStepSize_;
    nextStepSize_ = initialStepSize_;
  }

  stepSizeType_ = *stepSizeType;
  if (stepSizeType_ == STEP_TYPE_FIXED)
    currentStepSize_ = requestedStepSize_;
  else // STEP_TYPE_VARIABLE
    currentStepSize_ = nextStepSize_;

  // Limit the step size to the requested step size
  currentStepSize_ = std::min(requestedStepSize_, currentStepSize_);

  *stepSize = currentStepSize_;
  setStepControlState_(MID_STEP);
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::setCorrection(
     const StepperBase<Scalar>& stepper
    ,const RCP<const Thyra::VectorBase<Scalar> >& soln
    ,const RCP<const Thyra::VectorBase<Scalar> >& dx
    ,int solveStatus)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != MID_STEP, std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for SimpleStepControlStrategy<Scalar>::setCorrection()\n");
  x_ = soln;
  dx_ = dx;
  solveStatus_ = solveStatus;
  setStepControlState_(AFTER_CORRECTION);
}

template<class Scalar>
bool SimpleStepControlStrategy<Scalar>::acceptStep(
  const StepperBase<Scalar>& stepper, Scalar* value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for SimpleStepControlStrategy<Scalar>::completeStep()\n");

  if (solveStatus_ < 0 )
    return false;

  bool return_status = true;
  Scalar maxAbs_x  = std::max(Thyra::max(*x_),-Thyra::min(*x_));
  Scalar maxAbs_dx = std::max(Thyra::max(*dx_),-Thyra::min(*dx_));
  Scalar dx_tolerance = dxRelativeTolerance_ * maxAbs_x + dxAbsoluteTolerance_;
  if (maxAbs_dx > dx_tolerance) return_status = false;

  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if ( Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"acceptStep");
    *out << " max |*x_|    = " << maxAbs_x << "\n"
         << " dx_tolerance = " << dx_tolerance << "\n"
         << " max |*dx_|   = " << maxAbs_dx << "\n";
  }

  return(return_status);
}

template<class Scalar>
AttemptedStepStatusFlag SimpleStepControlStrategy<Scalar>::rejectStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for SimpleStepControlStrategy<Scalar>::completeStep()\n");

  setStepControlState_(READY_FOR_NEXT_STEP);

  using Teuchos::as;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"rejectStep");

  numStepFailures_ ++;
  if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) )
    *out << "numStepFailures_ = " << numStepFailures_ << "\n";
  if (numStepFailures_ > maxStepFailures_) {
    *out << "Rythmos_SimpleStepControlStrategy::rejectStep(...):  "
         << "Error: Too many step failures "
         << "(numStepFailures="<<numStepFailures_
         <<") > (maxStepFailures="<<maxStepFailures_<<")\n";
    return (REP_ERR_FAIL);
  }

  // Only update the time step if we are NOT running constant stepsize.
  if (stepSizeType_ == STEP_TYPE_VARIABLE) {
    nextStepSize_ *= stepSizeDecreaseFactor_;
    if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) ) {
      *out << "Rythmos_SimpleStepControl::rejectStep(...):  "
           << "  Step failure.  Reducing step size to "<< nextStepSize_ <<".\n";
    }
  } else {  // STEP_TYPE_FIXED
    if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) ) {
      *out << "Rythmos_SimpleStepControl::rejectStep(...):  "
           << "Error:  Step failure with fixed step size.\n";
    }
    return (REP_ERR_FAIL);
  }

  nextStepSize_ = std::max(nextStepSize_, minStepSize_);
  nextStepSize_ = std::min(nextStepSize_, maxStepSize_);

  AttemptedStepStatusFlag return_status = PREDICT_AGAIN;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "nextStepSize_ = " << nextStepSize_ << "\n";
  }

  return(return_status);
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::completeStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for SimpleStepControlStrategy<Scalar>::completeStep()\n");

  // Only update the time step if we are NOT running constant stepsize.
  if (stepSizeType_ == STEP_TYPE_VARIABLE) {
    // Only increase stepSize_ if no recent step failures.
    if (numStepFailures_ == 0) {
      nextStepSize_ *= stepSizeIncreaseFactor_;
    } else {
      // Keep nextStepSize_ constant until we have no recent step failures.
      nextStepSize_ = currentStepSize_;
      numStepFailures_ = std::max(numStepFailures_-1,0);
    }
  }
  nextStepSize_ = std::max(nextStepSize_, minStepSize_);
  nextStepSize_ = std::min(nextStepSize_, maxStepSize_);

  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if ( Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"completeStep_");
    *out << "nextStepSize_    = " << nextStepSize_ << "\n";
    *out << "numStepFailures_ = " << numStepFailures_ << "\n";
  }
  setStepControlState_(READY_FOR_NEXT_STEP);
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::as;

  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
       (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )    ) {
    out << this->description() << "::describe" << "\n";
  }
  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)) {
    out << "requestedStepSize_ = " << requestedStepSize_ << "\n";
    out << "currentStepSize_   = " << currentStepSize_ << "\n";
    out << "nextStepSize_      = " << nextStepSize_ << "\n";
    out << "stepSizeType_      = " << stepSizeType_ << "\n";
    out << "numStepFailures_   = " << numStepFailures_ << "\n";
  }
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList)
{
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_TEST_FOR_EXCEPT(paramList == Teuchos::null);
  paramList->validateParameters(*this->getValidParameters(),0);
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);

  initialStepSize_ = parameterList_->get(initialStepSizeName_,
                                         initialStepSizeDefault_);
  minStepSize_ = parameterList_->get(minStepSizeName_, minStepSizeDefault_);
  maxStepSize_ = parameterList_->get(maxStepSizeName_, maxStepSizeDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minStepSize_ <= maxStepSize_), std::logic_error,
    "Error:  (minStepSize="<<minStepSize_
    <<") > (maxStepSize="<<maxStepSize_<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    !((minStepSize_ <= initialStepSize_) && (initialStepSize_ <= maxStepSize_)),
    std::logic_error,
    "Error: Initial Step Size is not within min/max range.\n"
    << "  (minStepSize="<<minStepSize_
    <<") > (initialStepSize="<<initialStepSize_<<") or \n"
    << "  (maxStepSize="<<maxStepSize_
    <<") < (initialStepSize="<<initialStepSize_<<")\n");

  stepSizeIncreaseFactor_ = parameterList_->get(stepSizeIncreaseFactorName_,
                                                stepSizeIncreaseFactorDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(stepSizeIncreaseFactor_ > 0.0), std::logic_error,
    "Error:  (stepSizeIncreaseFactor="<<stepSizeIncreaseFactor_<<") <= 0.0\n");

  stepSizeDecreaseFactor_ = parameterList_->get(stepSizeDecreaseFactorName_,
                                                stepSizeDecreaseFactorDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(stepSizeDecreaseFactor_ > 0.0), std::logic_error,
    "Error:  (stepSizeDecreaseFactor="<<stepSizeDecreaseFactor_<<") <= 0.0\n");

  maxStepFailures_ = parameterList_->get(maxStepFailuresName_,
                                         maxStepFailuresDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(maxStepFailures_ >= 0), std::logic_error,
    "Error:  (maxStepFailures="<<maxStepFailures_<<") < 0\n");

  dxRelativeTolerance_ = parameterList_->get(dxRelativeToleranceName_,
                                             dxRelativeToleranceDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(dxRelativeTolerance_ >= 0.0), std::logic_error,
    "Error:  (dxRelativeTolerance="<<dxRelativeTolerance_<<") < 0.0\n");

  dxAbsoluteTolerance_ = parameterList_->get(dxAbsoluteToleranceName_,
                                             dxAbsoluteToleranceDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(dxAbsoluteTolerance_ >= 0.0), std::logic_error,
    "Error:  (dxAbsoluteTolerance="<<dxAbsoluteTolerance_<<") < 0.0\n");

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setParameterList");
  out->precision(15);

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "minStepSize_            = " << minStepSize_           << "\n";
    *out << "maxStepSize_            = " << maxStepSize_           << "\n";
    *out << "stepSizeIncreaseFactor_ = " << stepSizeIncreaseFactor_<< "\n";
    *out << "stepSizeDecreaseFactor_ = " << stepSizeDecreaseFactor_<< "\n";
    *out << "maxStepFailures_        = " << maxStepFailures_       << "\n";
    *out << "dxRelativeTolerance_    = " << dxRelativeTolerance_   << "\n";
    *out << "dxAbsoluteTolerance_    = " << dxAbsoluteTolerance_   << "\n";
  }
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
SimpleStepControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set(initialStepSizeName_,initialStepSizeDefault_, "Initial step size.");
    pl->set(minStepSizeName_, minStepSizeDefault_, "Minimum step size.");
    pl->set(maxStepSizeName_, maxStepSizeDefault_, "Maximum step size.");
    pl->set(stepSizeIncreaseFactorName_, stepSizeIncreaseFactorDefault_,
	    "Factor used to increase the step size after a successful step.");
    pl->set(stepSizeDecreaseFactorName_, stepSizeDecreaseFactorDefault_,
	    "Factor used to decrease the step size after a step failure.");
    pl->set(maxStepFailuresName_, maxStepFailuresDefault_,
      "The maximum number of step failures before exiting with an error.  "
      "The number of failure steps are carried between successful steps.");
    pl->set(dxRelativeToleranceName_, dxRelativeToleranceDefault_,
      "The allowable relative change in the solution for each step to "
      "pass.  The stepper solution status is also used to determine "
      "pass/fail.");
    pl->set(dxAbsoluteToleranceName_, dxAbsoluteToleranceDefault_,
      "The allowable absolute change in the solution for each step to "
      "pass.  The stepper solution status is also used to determine "
      "pass/fail.");

    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
SimpleStepControlStrategy<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
SimpleStepControlStrategy<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template<class Scalar>
void SimpleStepControlStrategy<Scalar>::setStepControlData(const StepperBase<Scalar>& stepper)
{
  if (stepControlState_ == UNINITIALIZED) initialize(stepper);
}

template<class Scalar>
bool SimpleStepControlStrategy<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
SimpleStepControlStrategy<Scalar>::cloneStepControlStrategyAlgorithm() const
{

  RCP<SimpleStepControlStrategy<Scalar> >
    stepControl = rcp(new SimpleStepControlStrategy<Scalar>());

  if (!is_null(parameterList_)) {
    stepControl->setParameterList(parameterList_);
  }

  return stepControl;
}

template<class Scalar>
int SimpleStepControlStrategy<Scalar>::getMaxOrder() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepControlState_ == UNINITIALIZED, std::logic_error,
      "Error, attempting to call getMaxOrder before initialization!\n"
      );
  return(maxOrder_);
}

//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_SIMPLE_STEP_CONTROL_STRATEGY_INSTANT(SCALAR) \
  template class SimpleStepControlStrategy< SCALAR >;


} // namespace Rythmos
#endif // Rythmos_SIMPLE_STEP_CONTROL_STRATEGY_DEF_H
