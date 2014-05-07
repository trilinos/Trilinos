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

#ifndef Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DEF_H
#define Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DEF_H

#include "Rythmos_FirstOrderErrorStepControlStrategy_decl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {

// Static members


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::initialStepSizeName_
= "Initial Step Size";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::initialStepSizeDefault_
= 1.0;


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::minStepSizeName_
= "Min Step Size";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::minStepSizeDefault_
= std::numeric_limits<Scalar>::min();


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::maxStepSizeName_
= "Max Step Size";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::maxStepSizeDefault_
= std::numeric_limits<Scalar>::max();


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::maxStepSizeIncreaseFactorName_
= "Max Step Size Increase Factor";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::maxStepSizeIncreaseFactorDefault_
= 1.5;


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::minStepSizeDecreaseFactorName_
= "Min Step Size Decrease Factor";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::minStepSizeDecreaseFactorDefault_
= 0.5;


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::maxStepFailuresName_
= "Maximum Number of Step Failures";

template<class Scalar>
const int
FirstOrderErrorStepControlStrategy<Scalar>::maxStepFailuresDefault_
= 100;


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::errorRelativeToleranceName_
= "Error Relative Tolerance";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::errorRelativeToleranceDefault_
= 1.0e-06;


template<class Scalar>
const std::string
FirstOrderErrorStepControlStrategy<Scalar>::errorAbsoluteToleranceName_
= "Error Absolute Tolerance";

template<class Scalar>
const double
FirstOrderErrorStepControlStrategy<Scalar>::errorAbsoluteToleranceDefault_
= 1.0e-12;


// Constructors

template<class Scalar>
void FirstOrderErrorStepControlStrategy<Scalar>::setStepControlState_(
  StepControlStrategyState newState)
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
StepControlStrategyState FirstOrderErrorStepControlStrategy<Scalar>::getCurrentState()
{
  return(stepControlState_);
}

template<class Scalar>
FirstOrderErrorStepControlStrategy<Scalar>::FirstOrderErrorStepControlStrategy()
  : stepControlState_(UNINITIALIZED),
    initialStepSize_(initialStepSizeDefault_),
    stepSizeType_(STEP_TYPE_VARIABLE),
    minStepSize_(minStepSizeDefault_),
    maxStepSize_(maxStepSizeDefault_),
    maxStepSizeIncreaseFactor_(maxStepSizeIncreaseFactorDefault_),
    minStepSizeDecreaseFactor_(minStepSizeDecreaseFactorDefault_),
    numStepFailures_(0),
    maxStepFailures_(maxStepFailuresDefault_),
    maxOrder_(1),
    errorRelativeTolerance_(errorRelativeToleranceDefault_),
    errorAbsoluteTolerance_(errorAbsoluteToleranceDefault_),
    solveStatus_(0)
{}

template<class Scalar>
void FirstOrderErrorStepControlStrategy<Scalar>::initialize(
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

  if (is_null(errWtVec_))
    errWtVec_ = createMember(stepper.get_x_space());
  setStepControlState_(BEFORE_FIRST_STEP);

  if (doTrace) {
    *out << "\nLeaving " << this->Teuchos::Describable::description()
         << "::initialize()...\n";
  }
}

template<class Scalar>
void FirstOrderErrorStepControlStrategy<Scalar>::setRequestedStepSize(
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
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::setRequestedStepSize()\n");

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
void FirstOrderErrorStepControlStrategy<Scalar>::nextStepSize(
  const StepperBase<Scalar>& stepper, Scalar* stepSize,
  StepSizeType* stepSizeType, int* order)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!((stepControlState_ == BEFORE_FIRST_STEP) ||
                               (stepControlState_ == MID_STEP) ||
                               (stepControlState_ == READY_FOR_NEXT_STEP) ),
                               std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::nextStepSize()\n");

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
void FirstOrderErrorStepControlStrategy<Scalar>::setCorrection(
     const StepperBase<Scalar>& stepper
    ,const RCP<const Thyra::VectorBase<Scalar> >& soln
    ,const RCP<const Thyra::VectorBase<Scalar> >& dx
    ,int solveStatus)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != MID_STEP, std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::setCorrection()\n");
  x_ = soln;
  dx_ = dx;
  solveStatus_ = solveStatus;
  setStepControlState_(AFTER_CORRECTION);
}

template<class Scalar>
bool FirstOrderErrorStepControlStrategy<Scalar>::acceptStep(
  const StepperBase<Scalar>& stepper, Scalar* value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::completeStep()\n");

  // Let's construct the weight here for now.  We will replace with it
  // with a errWtVecSet() in the future.
  Thyra::abs(*x_, Teuchos::ptrFromRef(*errWtVec_));
  Thyra::Vt_S(Teuchos::ptrFromRef(*errWtVec_), errorRelativeTolerance_);
  Thyra::Vp_S(Teuchos::ptrFromRef(*errWtVec_), errorAbsoluteTolerance_);
  reciprocal(*errWtVec_, Teuchos::ptrFromRef(*errWtVec_));
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // We square w because of how weighted norm_2 is computed.
  Vt_StV(Teuchos::ptrFromRef(*errWtVec_), ST::one(), *errWtVec_);
  // divide by N to get RMS norm
  int N = x_->space()->dim();
  Vt_S(Teuchos::ptrFromRef(*errWtVec_), Teuchos::as<Scalar>(1.0/N));
  double wrms = norm_2(*errWtVec_,(*dx_));
  stepSizeFactor_ = sqrt(2.0/wrms);     // Factor for 1st order. See
                                        // Gresho and Sani, "Incompressible
                                        // Flow and the Finite Element Method",
                                        // Vol. 1, 1998, p. 268.
  bool return_status = true;
  // If the error is too large (stepSizeFactor_ too small) OR
  // the solve failed, reduce the step size and try again.
  if ( (stepSizeFactor_ < minStepSizeDecreaseFactor_) || solveStatus_ < 0 )
    return_status = false;

  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  if ( Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ||
       (return_status == false &&
        Teuchos::as<int>(verbLevel) != Teuchos::as<int>(Teuchos::VERB_NONE)) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"acceptStep");
    *out << "\n wrms                       = " << wrms << "\n"
         << " stepSizeFactor_            = " << stepSizeFactor_ << "\n"
         << " solveStatus_               = " << solveStatus_ << "\n"
         << " minStepSizeDecreaseFactor_ = " << minStepSizeDecreaseFactor_
         << "\n";
  }

  return(return_status);
}

template<class Scalar>
AttemptedStepStatusFlag FirstOrderErrorStepControlStrategy<Scalar>::rejectStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::completeStep()\n");

  setStepControlState_(READY_FOR_NEXT_STEP);

  using Teuchos::as;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"rejectStep");

  numStepFailures_ ++;
  if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) )
    *out << "numStepFailures_ = " << numStepFailures_ << "\n";
  if (numStepFailures_ > maxStepFailures_) {
    *out << "Rythmos_FirstOrderErrorStepControlStrategy::rejectStep(...):  "
         << "Error: Too many step failures "
         << "(numStepFailures="<<numStepFailures_
         <<") > (maxStepFailures="<<maxStepFailures_<<")\n";
    return (REP_ERR_FAIL);
  }

  // Fail if we are running with fixed stepsize.
  if (stepSizeType_ == STEP_TYPE_FIXED) {
    if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) ) {
      *out << "Rythmos_FirstOrderErrorStepControl::rejectStep(...):  "
           << "Error:  Step failure with fixed step size.\n";
    }
    return (REP_ERR_FAIL);
  }

  nextStepSize_ = currentStepSize_*
                  std::max(stepSizeFactor_, minStepSizeDecreaseFactor_);
  nextStepSize_ = std::max(nextStepSize_, minStepSize_);
  nextStepSize_ = std::min(nextStepSize_, maxStepSize_);

  AttemptedStepStatusFlag return_status = PREDICT_AGAIN;

  if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) ) {
    *out << "Rythmos_FirstOrderErrorStepControl::rejectStep(...):  Step failure.\n"
         << "  Current step size is "<< currentStepSize_ <<".\n"
         << "  Reducing next step size to "<< nextStepSize_ <<".\n";
  }

  return(return_status);
}

template<class Scalar>
void FirstOrderErrorStepControlStrategy<Scalar>::completeStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepControlState_ != AFTER_CORRECTION,
     std::logic_error,
     "Error: Invalid state (stepControlState_=" << toString(stepControlState_)
     << ") for FirstOrderErrorStepControlStrategy<Scalar>::completeStep()\n");

  // Only update the time step if we are NOT running constant stepsize.
  if (stepSizeType_ == STEP_TYPE_VARIABLE) {
    // Only increase stepSize_ if no recent step failures.
    if (numStepFailures_ == 0) {
      nextStepSize_ *= std::max(minStepSizeDecreaseFactor_,
                         std::min(stepSizeFactor_, maxStepSizeIncreaseFactor_));
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
void FirstOrderErrorStepControlStrategy<Scalar>::describe(
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
void FirstOrderErrorStepControlStrategy<Scalar>::setParameterList(
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

  maxStepSizeIncreaseFactor_ =
    parameterList_->get(maxStepSizeIncreaseFactorName_,
                        maxStepSizeIncreaseFactorDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(maxStepSizeIncreaseFactor_ > 0.0), std::logic_error,
    "Error:  (maxStepSizeIncreaseFactor="
    <<maxStepSizeIncreaseFactor_<<") <= 0.0\n");

  minStepSizeDecreaseFactor_ =
    parameterList_->get(minStepSizeDecreaseFactorName_,
                        minStepSizeDecreaseFactorDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minStepSizeDecreaseFactor_ > 0.0), std::logic_error,
    "Error:  (minStepSizeDecreaseFactor="
    <<minStepSizeDecreaseFactor_<<") <= 0.0\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minStepSizeDecreaseFactor_ < maxStepSizeIncreaseFactor_),
    std::logic_error, "Error:  "
    "(minStepSizeDecreaseFactor="<<minStepSizeDecreaseFactor_<<") >="
    "(maxStepSizeIncreaseFactor="<<maxStepSizeIncreaseFactor_<<")\n");

  maxStepFailures_ = parameterList_->get(maxStepFailuresName_,
                                         maxStepFailuresDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(maxStepFailures_ >= 0), std::logic_error,
    "Error:  (maxStepFailures="<<maxStepFailures_<<") < 0\n");

  errorRelativeTolerance_ = parameterList_->get(errorRelativeToleranceName_,
                                                errorRelativeToleranceDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(errorRelativeTolerance_ >= 0.0), std::logic_error,
    "Error:  (errorRelativeTolerance="<<errorRelativeTolerance_<<") < 0.0\n");

  errorAbsoluteTolerance_ = parameterList_->get(errorAbsoluteToleranceName_,
                                                errorAbsoluteToleranceDefault_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(errorAbsoluteTolerance_ >= 0.0), std::logic_error,
    "Error:  (errorAbsoluteTolerance="<<errorAbsoluteTolerance_<<") < 0.0\n");

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setParameterList");
  out->precision(15);

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "minStepSize_               = " << minStepSize_              <<"\n";
    *out << "maxStepSize_               = " << maxStepSize_              <<"\n";
    *out << "maxStepSizeIncreaseFactor_ = " << maxStepSizeIncreaseFactor_<<"\n";
    *out << "minStepSizeDecreaseFactor_ = " << minStepSizeDecreaseFactor_<<"\n";
    *out << "maxStepFailures_           = " << maxStepFailures_          <<"\n";
    *out << "errorRelativeTolerance_    = " << errorRelativeTolerance_   <<"\n";
    *out << "errorAbsoluteTolerance_    = " << errorAbsoluteTolerance_   <<"\n";
  }
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
FirstOrderErrorStepControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set(initialStepSizeName_,initialStepSizeDefault_, "Initial step size.");
    pl->set(minStepSizeName_, minStepSizeDefault_, "Minimum step size.");
    pl->set(maxStepSizeName_, maxStepSizeDefault_, "Maximum step size.");
    pl->set(maxStepSizeIncreaseFactorName_, maxStepSizeIncreaseFactorDefault_,
	    "The maximum factor to increase the step size after a successful step.");
    pl->set(minStepSizeDecreaseFactorName_, minStepSizeDecreaseFactorDefault_,
	    "The minimum allowable factor to decrease the step size.  If the "
      "stepSizeFactor_ is below this, the current step is considered a "
      "failure and retried with `" + minStepSizeDecreaseFactorName_ +
      "' applied.");
    pl->set(maxStepFailuresName_, maxStepFailuresDefault_,
      "The maximum number of step failures before exiting with an error.  "
      "The number of failure steps are carried between successful steps.");
    pl->set(errorRelativeToleranceName_, errorRelativeToleranceDefault_,
      "The allowable relative change in the error (the difference between "
      "the solution at the end of the step and the predicted solution) for "
      "each step to pass.  The stepper solution status is also used to "
      "determine pass/fail.");
    pl->set(errorAbsoluteToleranceName_, errorAbsoluteToleranceDefault_,
      "The allowable absolute change in the error (the difference between "
      "the solution at the end of the step and the predicted solution) for "
      "each step to pass.  The stepper solution status is also used to "
      "determine pass/fail.");

    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
FirstOrderErrorStepControlStrategy<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
FirstOrderErrorStepControlStrategy<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template<class Scalar>
void FirstOrderErrorStepControlStrategy<Scalar>::setStepControlData(
  const StepperBase<Scalar>& stepper)
{
  if (stepControlState_ == UNINITIALIZED) initialize(stepper);
}

template<class Scalar>
bool FirstOrderErrorStepControlStrategy<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
FirstOrderErrorStepControlStrategy<Scalar>::cloneStepControlStrategyAlgorithm() const
{
  RCP<FirstOrderErrorStepControlStrategy<Scalar> >
    stepControl = rcp(new FirstOrderErrorStepControlStrategy<Scalar>());

  if (!is_null(parameterList_)) {
    stepControl->setParameterList(parameterList_);
  }

  return stepControl;
}

template<class Scalar>
int FirstOrderErrorStepControlStrategy<Scalar>::getMaxOrder() const
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

#define RYTHMOS_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_INSTANT(SCALAR) \
  template class FirstOrderErrorStepControlStrategy< SCALAR >;


} // namespace Rythmos
#endif // Rythmos_FIRSTORDERERROR_STEP_CONTROL_STRATEGY_DEF_H
