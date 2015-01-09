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

#ifndef Rythmos_IMPLICITBDF_STEPPER_RAMPING_STEP_CONTROL_DEF_H
#define Rythmos_IMPLICITBDF_STEPPER_RAMPING_STEP_CONTROL_DEF_H

#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Rythmos {

template<class Scalar>
ImplicitBDFStepperRampingStepControl<Scalar>::
ImplicitBDFStepperRampingStepControl() :
  stepControlState_(UNINITIALIZED)
{

}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setStepControlState_(
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
StepControlStrategyState
ImplicitBDFStepperRampingStepControl<Scalar>::getCurrentState()
{
  return(stepControlState_);
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::updateCoeffs_()
{
  TEUCHOS_TEST_FOR_EXCEPT(!((stepControlState_ == BEFORE_FIRST_STEP) ||
                            (stepControlState_ == READY_FOR_NEXT_STEP)));

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "updateCoeffs_() is not implemented!");
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::initialize(
  const StepperBase<Scalar>& stepper)
{
  // Initialize can be called from the stepper when setInitialCondition
  // is called.
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  // Set initial time:
  TimeRange<Scalar> stepperRange = stepper.getTimeRange();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !stepperRange.isValid(),
      std::logic_error,
      "Error, Stepper does not have valid time range for initialization "
      "of ImplicitBDFStepperRampingStepControl!\n");

  if (is_null(parameterList_)) {
    RCP<Teuchos::ParameterList> emptyParameterList =
      Teuchos::rcp(new Teuchos::ParameterList);
    this->setParameterList(emptyParameterList);
  }

  if (is_null(errWtVecCalc_)) {
    RCP<ImplicitBDFStepperErrWtVecCalc<Scalar> > IBDFErrWtVecCalc =
      rcp(new ImplicitBDFStepperErrWtVecCalc<Scalar>());
    errWtVecCalc_ = IBDFErrWtVecCalc;
  }

  stepControlState_ = UNINITIALIZED;

  requestedStepSize_ = Scalar(-1.0);
  currentStepSize_ = initialStepSize_;
  currentOrder_ = 1;
  nextStepSize_ = initialStepSize_;
  nextOrder_ = 1;
  numberOfSteps_ = 0;
  totalNumberOfFailedSteps_ = 0;
  countOfConstantStepsAfterFailure_ = 0;

  if (is_null(delta_)) {
    delta_ = createMember(stepper.get_x_space());
  }
  if (is_null(errWtVec_)) {
    errWtVec_ = createMember(stepper.get_x_space());
  }
  V_S(delta_.ptr(),ST::zero());

  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"initialize");
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
  }

  setStepControlState_(BEFORE_FIRST_STEP);

}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setRequestedStepSize(
    const StepperBase<Scalar>& stepper,
    const Scalar& stepSize,
    const StepSizeType& stepSizeType)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_TEST_FOR_EXCEPT(!((stepControlState_ == UNINITIALIZED) ||
                            (stepControlState_ == BEFORE_FIRST_STEP) ||
                            (stepControlState_ == READY_FOR_NEXT_STEP) ||
                            (stepControlState_ == MID_STEP)));

  TEUCHOS_TEST_FOR_EXCEPTION(
      ((stepSizeType == STEP_TYPE_FIXED) && (stepSize == ST::zero())),
      std::logic_error,
      "Error, step size type == STEP_TYPE_FIXED, "
      "but requested step size == 0!\n");

  bool didInitialization = false;
  if (stepControlState_ == UNINITIALIZED) {
    initialize(stepper);
    didInitialization = true;
  }

  // errWtVecSet_ is called during initialize
  if (!didInitialization) {
    const ImplicitBDFStepper<Scalar>& implicitBDFStepper =
      Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
    const Thyra::VectorBase<Scalar>& xHistory =
      implicitBDFStepper.getxHistory(0);
    errWtVecCalc_->errWtVecSet(&*errWtVec_,xHistory,relErrTol_,absErrTol_);
  }

  requestedStepSize_ = stepSize;
  stepSizeType_ = stepSizeType;
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::nextStepSize(
  const StepperBase<Scalar>& stepper, Scalar* stepSize,
  StepSizeType* stepSizeType, int* order)
{
  TEUCHOS_TEST_FOR_EXCEPT(!((stepControlState_ == BEFORE_FIRST_STEP) ||
         (stepControlState_ == MID_STEP) ||
         (stepControlState_ == READY_FOR_NEXT_STEP) )
        );

  if (stepControlState_ == BEFORE_FIRST_STEP) {
    nextStepSize_ = initialStepSize_;
    nextOrder_ = 1;
  }

  // Now starting a step - rotate next values into current values
  if (stepSizeType_ == STEP_TYPE_FIXED)
    currentStepSize_ = requestedStepSize_;
  else
    currentStepSize_ = nextStepSize_;

  currentOrder_ = nextOrder_;

  // Limit the step size to the requested step size
  currentStepSize_ = std::min(requestedStepSize_, currentStepSize_);

  *stepSize = currentStepSize_;
  *stepSizeType = stepSizeType_;
  *order = currentOrder_;

  if (stepControlState_ != MID_STEP) {
    setStepControlState_(MID_STEP);
  }

  // Output
  if (doOutput_(Teuchos::VERB_MEDIUM)){
    Teuchos::FancyOStream& out = *this->getOStream();
    Teuchos::OSTab ostab1(out,2,"** nextStepSize_ **");
    out << "Values returned to stepper:" << std::endl;
    Teuchos::OSTab ostab2(out,2,"** nextStepSize_ **");
    out << "currentStepSize_ = " << currentStepSize_ << std::endl;
    out << "currentOrder_ = " << currentOrder_ << std::endl;
    out << "requestedStepSize_ = " << requestedStepSize_ << std::endl;
  }

}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setCorrection(
     const StepperBase<Scalar>& stepper
    ,const RCP<const Thyra::VectorBase<Scalar> >& soln
    ,const RCP<const Thyra::VectorBase<Scalar> >& ee
    ,int solveStatus)
{
  TEUCHOS_TEST_FOR_EXCEPT(stepControlState_ != MID_STEP);

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(ee), std::logic_error,
    "Error, ee == Teuchos::null!\n");

  ee_ = ee;

  newtonConvergenceStatus_ = solveStatus;

  if ( doOutput_(Teuchos::VERB_MEDIUM)  && newtonConvergenceStatus_ < 0) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"setCorrection");
    *out << "\nImplicitBDFStepperRampingStepControl::setCorrection(): "
         << "Nonlinear Solver Failed!\n";
  }

  setStepControlState_(AFTER_CORRECTION);
}

template<class Scalar>
bool ImplicitBDFStepperRampingStepControl<Scalar>::acceptStep(
  const StepperBase<Scalar>& stepper, Scalar* LETValue)
{
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);


  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::OSTab ostab(out,1,"acceptStep");
    *out << "ee_ = " << std::endl;
    ee_->describe(*out,verbLevel);
    *out << "errWtVec_ = " << std::endl;
    errWtVec_->describe(*out,verbLevel);
  }

  Scalar enorm = wRMSNorm_(*errWtVec_,*ee_);

  Scalar LET = ck_ * enorm;

  if (LETValue) {
    *LETValue = LET;
    *LETValue = Scalar(0.0);
  }

  if (newtonConvergenceStatus_ < 0)
    return false;

  bool return_status = false;

  if (LET < ST::one() || !useLETToDetermineConvergence_)
    return_status = true;

  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"acceptStep");
    *out << "return_status = " << return_status << std::endl;
    *out << "Local Truncation Error Check: (ck*enorm) < 1:  (" << LET
         << ") <?= 1" << std::endl;
    if ( doOutput_(Teuchos::VERB_EXTREME) ) {
      *out << "ck_ = " << ck_ << std::endl;
      *out << "enorm = " << enorm << std::endl;
    }
  }

  return(return_status);
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::completeStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();

    Teuchos::OSTab ostab1(out,2,"completeStep_");
    *out << "\n** Begin completeStep() **" << std::endl;
    Teuchos::OSTab ostab2(out,2,"** Begin completeStep_ **");
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "numConstantSteps_ = " << numConstantSteps_ << std::endl;
    *out << "currentStepSize_ = " << currentStepSize_ << std::endl;
    *out << "nextStepSize_ = " << nextStepSize_ << std::endl;
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "nextOrder_ = " << nextOrder_ << std::endl;
    *out << "stepSizeIncreaseFactor_ = " << stepSizeIncreaseFactor_ <<std::endl;
    *out << "countOfConstantStepsAfterFailure_ = "
         << countOfConstantStepsAfterFailure_ << std::endl;
  }

  numberOfSteps_ ++;

  if (countOfConstantStepsAfterFailure_ > 0) {
    // We track the number of consecutive time step failures so that
    // if we have a bunch of nonlinear failures, lets keep the time
    // step constant for a while before we start to ramp again.  This
    // keeps us from oscillating between ramping and cutting step
    // sizes and wasting resources.

    nextStepSize_ = currentStepSize_;
    nextOrder_ = currentOrder_;

    // Decrement failure counter
    countOfConstantStepsAfterFailure_ =
      std::max( (countOfConstantStepsAfterFailure_ - 1), 0);

    if ( doOutput_(Teuchos::VERB_HIGH) ) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"completeStep_");
      *out << "\nNext Step Size held constant due to previous failed steps!\n";
      *out << "countOfConstantStepsAfterFailure_ = "
           << countOfConstantStepsAfterFailure_ << std::endl;
    }
  }
  else {

    // Phase 1: Constant step size at 1st order
    if (numberOfSteps_ < numConstantSteps_) {
      if (currentStepSize_ <  initialStepSize_)
        nextStepSize_ = std::min(initialStepSize_,
                                 currentStepSize_ * stepSizeIncreaseFactor_);
      nextOrder_ = 1;
    }
    // Phase 2: Constant step size, ramping the order
    else if (currentOrder_ < maxOrder_) {
      if (currentStepSize_ <  initialStepSize_)
        nextStepSize_ = std::min(initialStepSize_,
                                 currentStepSize_ * stepSizeIncreaseFactor_);
      else
        nextStepSize_ = currentStepSize_;

      nextOrder_ = currentOrder_ + 1;
    }
    // Phase 3: Ramping dt to max step size, highest order
    else if ( (numberOfSteps_ >= numConstantSteps_) &&
              (currentOrder_ == maxOrder_)             ) {
      nextStepSize_ = std::min(maxStepSize_,
                               currentStepSize_ * stepSizeIncreaseFactor_);
      nextOrder_ = maxOrder_;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "RampingStepControlStrategy logic is broken. Please contact "
        "developers. Aborting run!");
    }

    if (restrictStepSizeByNumberOfNonlinearIterations_) {
      const Rythmos::ImplicitBDFStepper<Scalar>* ibdfStepper =
	dynamic_cast<const Rythmos::ImplicitBDFStepper<Scalar>* >(&stepper);
      TEUCHOS_ASSERT(ibdfStepper != NULL);
      TEUCHOS_ASSERT(nonnull(ibdfStepper->getNonlinearSolveStatus().extraParameters));
      int numberOfNonlinearIterations = ibdfStepper->getNonlinearSolveStatus().extraParameters->template get<int>("Number of Iterations");
      if (numberOfNonlinearIterations >= numberOfNonlinearIterationsForStepSizeRestriction_) {
	nextStepSize_ = currentStepSize_;
      }
    }


  } // if (countOfConstantStepsAfterFailure_ > 0)

  setStepControlState_(READY_FOR_NEXT_STEP);

  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab1(out,2,"** completeStep_ **");
    *out << "** End of completeStep() **" << std::endl;
    Teuchos::OSTab ostab2(out,2,"** End completeStep_ **");
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "numConstantSteps_ = " << numConstantSteps_ << std::endl;
    *out << "currentStepSize_ = " << currentStepSize_ << std::endl;
    *out << "nextStepSize_ = " << nextStepSize_ << std::endl;
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "nextOrder_ = " << nextOrder_ << std::endl;
    *out << "stepSizeIncreaseFactor_ = " << stepSizeIncreaseFactor_ <<std::endl;
    *out << "countOfConstantStepsAfterFailure_ = "
         << countOfConstantStepsAfterFailure_ << std::endl;
  }
}

template<class Scalar>
AttemptedStepStatusFlag
ImplicitBDFStepperRampingStepControl<Scalar>::rejectStep(
  const StepperBase<Scalar>& stepper)
{
  TEUCHOS_TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);

  using Teuchos::as;

  ++totalNumberOfFailedSteps_;
  ++countOfConstantStepsAfterFailure_;

  // If time step size is already at the min time step, then quit
  if (currentStepSize_ <= minStepSize_)
    return (REP_ERR_FAIL);

  // Otherwise, cut the time step and keep order the same
  nextStepSize_ = std::max(minStepSize_,
                           (currentStepSize_ * stepSizeDecreaseFactor_) );

  setStepControlState_(READY_FOR_NEXT_STEP);

  return(PREDICT_AGAIN);
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::as;

  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
    (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )
    ) {
    out << this->description() << "::describe" << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)) {
    out << "currentStepSize_ = " << currentStepSize_ << std::endl;
    out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM)) {
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
    out << "ee_ = ";
    if (ee_ == Teuchos::null) {
      out << "Teuchos::null" << std::endl;
    } else {
      ee_->describe(out,verbLevel);
    }
    out << "delta_ = ";
    if (delta_ == Teuchos::null) {
      out << "Teuchos::null" << std::endl;
    } else {
      delta_->describe(out,verbLevel);
    }
    out << "errWtVec_ = ";
    if (errWtVec_ == Teuchos::null) {
      out << "Teuchos::null" << std::endl;
    } else {
      errWtVec_->describe(out,verbLevel);
    }
  }
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_TEST_FOR_EXCEPT(paramList == Teuchos::null);

  parameterList_ = Teuchos::parameterList(*paramList);

  parameterList_->validateParametersAndSetDefaults(*this->getValidParameters());

  Teuchos::ParameterList& p = *parameterList_;

  numConstantSteps_ = p.get<int>("Number of Constant First Order Steps");
  initialStepSize_ = p.get<Scalar>("Initial Step Size");
  minStepSize_ = p.get<Scalar>("Min Step Size");
  maxStepSize_ = p.get<Scalar>("Max Step Size");
  stepSizeIncreaseFactor_ = p.get<Scalar>("Step Size Increase Factor");
  stepSizeDecreaseFactor_ = p.get<Scalar>("Step Size Decrease Factor");

  minOrder_ = p.get<int>("Min Order");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !((1 <= minOrder_) && (minOrder_ <= 5)), std::logic_error,
      "Error, minOrder_ = " << minOrder_ << " is not in range [1,5]!\n"
      );
  maxOrder_ = p.get<int>("Max Order");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !((1 <= maxOrder_) && (maxOrder_ <= 5)), std::logic_error,
      "Error, maxOrder_ = " << maxOrder_ << " is not in range [1,5]!\n"
      );

  absErrTol_ = p.get<Scalar>("Absolute Error Tolerance");
  relErrTol_ = p.get<Scalar>("Relative Error Tolerance");

  {
    std::string let_acceptance =
      p.get<std::string>("Use LET To Determine Step Acceptance");
    useLETToDetermineConvergence_ = (let_acceptance == "TRUE");

    // Currently the using LET for step acceptance is not supported
    // since we can't calculate the LETValue. Once this is
    // implemented, delete the line below.
    TEUCHOS_TEST_FOR_EXCEPTION(useLETToDetermineConvergence_, std::logic_error,
      "Error - the flag \"Use LET To Determine Step Acceptance\" is set to "
      "\"TRUE\" but the local error computation is currently not supported.  "
      "Please set this flag to \"FALSE\" for now.");
  }

  if (p.get<std::string>("Restrict Step Size Increase by Number of Nonlinear Iterations") == "TRUE")
    restrictStepSizeByNumberOfNonlinearIterations_ = true;
  else if (p.get<std::string>("Restrict Step Size Increase by Number of Nonlinear Iterations") == "FALSE") 
    restrictStepSizeByNumberOfNonlinearIterations_ = false;
    
  numberOfNonlinearIterationsForStepSizeRestriction_ = 
    p.get<int>("Number of Nonlinear Iterations for Step Size Restriction");

  if ( doOutput_(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"setParameterList");
    out->precision(15);
    *out << "minOrder_ = " << minOrder_ << std::endl;
    *out << "maxOrder_ = " << maxOrder_ << std::endl;
    *out << "relErrTol_  = " << relErrTol_  << std::endl;
    *out << "absErrTol_  = " << absErrTol_  << std::endl;
    *out << "stepSizeType = " << stepSizeType_  << std::endl;
    *out << "stopTime_  = " << stopTime_  << std::endl;
  }

}

template<class Scalar>
RCP<const Teuchos::ParameterList>
ImplicitBDFStepperRampingStepControl<Scalar>::getValidParameters() const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  static RCP<ParameterList> p;

  if (is_null(p)) {

    p = rcp(new ParameterList);

    p->set<int>("Number of Constant First Order Steps", 10,
      "Number of constant steps to take before handing control to "
      "variable stepper.");
    p->set<Scalar>("Initial Step Size", Scalar(1.0e-3),
      "Initial time step size and target step size to take during the "
      "initial constant step phase (could be reduced due to step failures).");
    p->set<Scalar>("Min Step Size", Scalar(1.0e-7), "Minimum time step size.");
    p->set<Scalar>("Max Step Size", Scalar(1.0), "Maximum time step size.");
    p->set<Scalar>("Step Size Increase Factor", Scalar(1.2),
      "Time step growth factor used after a successful time step. dt_{n+1} = "
      "(increase factor) * dt_n");
    p->set<Scalar>("Step Size Decrease Factor", Scalar(0.5),
      "Time step reduction factor used for a failed time step. dt_{n+1} = "
      "(decrease factor) * dt_n");
    p->set<int>("Min Order", 1, "Minimum order to run at.");
    p->set<int>("Max Order", 5, "Maximum order to run at.");
    p->set<Scalar>("Absolute Error Tolerance", Scalar(1.0e-5),
      "abstol value used in WRMS calculation.");
    p->set<Scalar>("Relative Error Tolerance", Scalar(1.0e-3),
      "reltol value used in WRMS calculation.");
    Teuchos::setStringToIntegralParameter<int>(
      "Use LET To Determine Step Acceptance",
      "FALSE",
      "If set to TRUE, then acceptance of step dependes on LET in addition "
      "to Nonlinear solver converging.",
      Teuchos::tuple<std::string>("TRUE","FALSE"),
      p.get());
    Teuchos::setStringToIntegralParameter<int>(
      "Restrict Step Size Increase by Number of Nonlinear Iterations",
      "FALSE",
      "If set to TRUE, then the step size will not be allowed to increase "
      "if the number of nonlinear iterations was greater than or equal to the "
      "specified value.",
      Teuchos::tuple<std::string>("TRUE","FALSE"),
      p.get());
    p->set<int>("Number of Nonlinear Iterations for Step Size Restriction",
		2,
		"If \" Restrct Step Size Increase by Number of Nonlinear Iterations\" is "
                "true, the step size will not be allowed to increase if the number of nonlinear "
		"iterations was greater than or equal to the specified value.");
  }

  return (p);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperRampingStepControl<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperRampingStepControl<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setStepControlData(
  const StepperBase<Scalar>& stepper)
{
  if (stepControlState_ == UNINITIALIZED) {
    initialize(stepper);
  }
  const ImplicitBDFStepper<Scalar>& bdfstepper =
    Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
  int desiredOrder = bdfstepper.getOrder();
  TEUCHOS_TEST_FOR_EXCEPT(!((1 <= desiredOrder) &&
                            (desiredOrder <= maxOrder_)));
  if (stepControlState_ == BEFORE_FIRST_STEP) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        desiredOrder > 1,
        std::logic_error,
        "Error, this ImplicitBDF stepper has not taken a step yet, so it "
        "cannot take a step of order " << desiredOrder << " > 1!\n");
  }
  TEUCHOS_TEST_FOR_EXCEPT(!(desiredOrder <= currentOrder_+1));
  currentOrder_ = desiredOrder;

  if ( doOutput_(Teuchos::VERB_EXTREME) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"setStepControlData");
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
}

template<class Scalar>
bool ImplicitBDFStepperRampingStepControl<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
ImplicitBDFStepperRampingStepControl<Scalar>::cloneStepControlStrategyAlgorithm() const
{

  RCP<ImplicitBDFStepperRampingStepControl<Scalar> > stepControl =
    rcp(new ImplicitBDFStepperRampingStepControl<Scalar>());

  if (!is_null(parameterList_)) {
    stepControl->setParameterList(parameterList_);
  }

  return stepControl;
}

template<class Scalar>
void ImplicitBDFStepperRampingStepControl<Scalar>::setErrWtVecCalc(
  const RCP<ErrWtVecCalcBase<Scalar> >& errWtVecCalc)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(errWtVecCalc));
  errWtVecCalc_ = errWtVecCalc;
}

template<class Scalar>
RCP<const ErrWtVecCalcBase<Scalar> >
ImplicitBDFStepperRampingStepControl<Scalar>::getErrWtVecCalc() const
{
  return(errWtVecCalc_);
}

template<class Scalar>
Scalar ImplicitBDFStepperRampingStepControl<Scalar>::wRMSNorm_(
    const Thyra::VectorBase<Scalar>& weight,
    const Thyra::VectorBase<Scalar>& vector) const
{
  return(norm_2(weight,vector));
}

template<class Scalar>
int ImplicitBDFStepperRampingStepControl<Scalar>::getMinOrder() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepControlState_ == UNINITIALIZED, std::logic_error,
      "Error, attempting to call getMinOrder before intiialization!\n"
      );
  return(minOrder_);
}

template<class Scalar>
int ImplicitBDFStepperRampingStepControl<Scalar>::getMaxOrder() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepControlState_ == UNINITIALIZED, std::logic_error,
      "Error, attempting to call getMaxOrder before initialization!\n"
      );
  return(maxOrder_);
}

template<class Scalar>
bool ImplicitBDFStepperRampingStepControl<Scalar>::doOutput_(
  Teuchos::EVerbosityLevel verbLevel)
{
  Teuchos::EVerbosityLevel currentObjectVerbLevel = this->getVerbLevel();

  if ( Teuchos::as<int>(currentObjectVerbLevel) >= Teuchos::as<int>(verbLevel) )
    return true;

  return false;
}

template<class Scalar>
int ImplicitBDFStepperRampingStepControl<Scalar>::numberOfSteps() const
{
  return numberOfSteps_;
}

template<class Scalar>
int ImplicitBDFStepperRampingStepControl<Scalar>::numberOfFailedSteps() const
{
  return totalNumberOfFailedSteps_;
}

template<class Scalar>
Scalar ImplicitBDFStepperRampingStepControl<Scalar>::currentStepSize() const
{
  return currentStepSize_;
}

template<class Scalar>
int ImplicitBDFStepperRampingStepControl<Scalar>::currentOrder() const
{
  return currentOrder_;
}


//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_IMPLICITBDF_STEPPER_RAMPING_STEPCONTROL_INSTANT(SCALAR) \
  template class ImplicitBDFStepperRampingStepControl< SCALAR >;


} // namespace Rythmos
#endif // Rythmos_IMPLICITBDF_STEPPER_RAMPING_STEP_CONTROL_DEF_H

