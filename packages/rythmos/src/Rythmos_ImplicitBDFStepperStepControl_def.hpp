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

#ifndef Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DEF_H
#define Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DEF_H

#include "Rythmos_ImplicitBDFStepperStepControl_decl.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ImplicitBDFStepperErrWtVecCalc.hpp"

namespace Rythmos {

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setStepControlState_(StepControlStrategyState newState)
{
  if (stepControlState_ == UNINITIALIZED) {
    TEST_FOR_EXCEPT(newState != BEFORE_FIRST_STEP);
  } else if (stepControlState_ == BEFORE_FIRST_STEP) {
    TEST_FOR_EXCEPT(newState != MID_STEP);
  } else if (stepControlState_ == MID_STEP) {
    TEST_FOR_EXCEPT(newState != AFTER_CORRECTION);
  } else if (stepControlState_ == AFTER_CORRECTION) {
    TEST_FOR_EXCEPT(newState != READY_FOR_NEXT_STEP);
    checkReduceOrderCalled_ = false;
  } else if (stepControlState_ == READY_FOR_NEXT_STEP) {
    TEST_FOR_EXCEPT(newState != MID_STEP);
  }
  stepControlState_ = newState;
}

template<class Scalar>
StepControlStrategyState ImplicitBDFStepperStepControl<Scalar>::getCurrentState()
{
  return(stepControlState_);
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::updateCoeffs_() 
{
  TEST_FOR_EXCEPT(!((stepControlState_ == BEFORE_FIRST_STEP) || (stepControlState_ == READY_FOR_NEXT_STEP)));
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((hh_ != usedStep_) || (currentOrder_ != usedOrder_)) {
    nscsco_ = 0;
  }
  nscsco_ = std::min(nscsco_+1,usedOrder_+2);
  if (currentOrder_+1 >= nscsco_) {
    alpha_[0] = ST::one();
    Scalar temp1 = hh_;
    sigma_[0] = ST::one();
    gamma_[0] = ST::zero();
    for (int i=1;i<=currentOrder_;++i) {
      Scalar temp2 = psi_[i-1];
      psi_[i-1] = temp1;
      temp1 = temp2 + hh_;
      alpha_[i] = hh_/temp1;
      sigma_[i] = Scalar(i+1)*sigma_[i-1]*alpha_[i];
      gamma_[i] = gamma_[i-1]+alpha_[i-1]/hh_;
    }
    psi_[currentOrder_] = temp1;
    alpha_s_ = ST::zero();
    alpha_0_ = ST::zero();
    for (int i=0;i<currentOrder_;++i) {
      alpha_s_ = alpha_s_ - Scalar(ST::one()/(i+ST::one()));
      alpha_0_ = alpha_0_ - alpha_[i];
    }
    cj_ = -alpha_s_/hh_;
    ck_ = std::abs(alpha_[currentOrder_]+alpha_s_-alpha_0_);
    ck_ = std::max(ck_,alpha_[currentOrder_]);
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"updateCoeffs_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    for (int i=0;i<=maxOrder_;++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "sigma_[" << i << "] = " << sigma_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
      *out << "alpha_s_ = " << alpha_s_ << std::endl;
      *out << "alpha_0_ = " << alpha_0_ << std::endl;
      *out << "cj_ = " << cj_ << std::endl;
      *out << "ck_ = " << ck_ << std::endl;
    }
  }
}

template<class Scalar>
ImplicitBDFStepperStepControl<Scalar>::ImplicitBDFStepperStepControl()
{
  defaultInitializeAllData_();
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::initialize(const StepperBase<Scalar>& stepper)
{
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool doTrace = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH));
  Teuchos::OSTab ostab(out,1,"initialize");

  if (doTrace) {
    *out
      << "\nEntering " << this->Teuchos::Describable::description()
      << "::initialize()...\n";
  }

  // Set initial time:
  TimeRange<Scalar> stepperRange = stepper.getTimeRange();
  TEST_FOR_EXCEPTION(
      !stepperRange.isValid(),
      std::logic_error,
      "Error, Stepper does not have valid time range for initialization of ImplicitBDFStepperStepControl!\n"
      );
  time_ = stepperRange.upper();

  if (parameterList_ == Teuchos::null) {
    RCP<Teuchos::ParameterList> emptyParameterList = Teuchos::rcp(new Teuchos::ParameterList);
    this->setParameterList(emptyParameterList);
  }

  if (is_null(errWtVecCalc_)) {
    RCP<ImplicitBDFStepperErrWtVecCalc<Scalar> > IBDFErrWtVecCalc = rcp(new ImplicitBDFStepperErrWtVecCalc<Scalar>());
    errWtVecCalc_ = IBDFErrWtVecCalc;
  }

  // 08/22/07 initialize can be called from the stepper when setInitialCondition is called.
  checkReduceOrderCalled_ = false;
  stepControlState_ = UNINITIALIZED;

  currentOrder_=1; // Current order of integration
  oldOrder_=1; // previous order of integration
  usedOrder_=1;  // order used in current step (used after currentOrder_ is updated)
  alpha_s_=Scalar(-ST::one());  // $\alpha_s$ fixed-leading coefficient of this BDF method
  alpha_.clear();  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
  // note:   $h_n$ = current step size, n = current time step
  gamma_.clear();  // calculate time derivative of history array for predictor 
  psi_.clear();    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
  sigma_.clear();  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  Scalar zero = ST::zero();
  for (int i=0 ; i<=maxOrder_ ; ++i) {
    alpha_.push_back(zero);
    gamma_.push_back(zero);
    psi_.push_back(zero);
    sigma_.push_back(zero);
  }
  alpha_0_=zero;   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  cj_=zero;      // $-\alpha_s/h_n$ coefficient used in local error test
  ck_=zero;      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  hh_=zero;
  numberOfSteps_=0;   // number of total time integration steps taken
  nef_=0;
  usedStep_=zero;
  nscsco_=0;
  Ek_=zero;
  Ekm1_=zero;
  Ekm2_=zero;
  Ekp1_=zero;
  Est_=zero;
  Tk_=zero;
  Tkm1_=zero;
  Tkm2_=zero;
  Tkp1_=zero;
  newOrder_=currentOrder_;
  initialPhase_=true;

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "oldOrder_ = " << oldOrder_ << std::endl;
    *out << "usedOrder_ = " << usedOrder_ << std::endl;
    *out << "alpha_s_ = " << alpha_s_ << std::endl;
    for (int i=0 ; i<=maxOrder_ ; ++i) {
      *out << "alpha_[" << i << "] = " << alpha_[i] << std::endl;
      *out << "gamma_[" << i << "] = " << gamma_[i] << std::endl;
      *out << "psi_[" << i << "] = " << psi_[i] << std::endl;
      *out << "sigma_[" << i << "] = " << sigma_[i] << std::endl;
    }
    *out << "alpha_0_ = " << alpha_0_ << std::endl;
    *out << "cj_ = " << cj_ << std::endl;
    *out << "ck_ = " << ck_ << std::endl;
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
    *out << "usedStep_ = " << usedStep_ << std::endl;
    *out << "nscsco_ = " << nscsco_ << std::endl;
    *out << "Ek_ = " << Ek_ << std::endl;
    *out << "Ekm1_ = " << Ekm1_ << std::endl;
    *out << "Ekm2_ = " << Ekm2_ << std::endl;
    *out << "Ekp1_ = " << Ekp1_ << std::endl;
    *out << "Est_ = " << Est_ << std::endl;
    *out << "Tk_ = " << Tk_ << std::endl;
    *out << "Tkm1_ = " << Tkm1_ << std::endl;
    *out << "Tkm2_ = " << Tkm2_ << std::endl;
    *out << "Tkp1_ = " << Tkp1_ << std::endl;
    *out << "newOrder_ = " << newOrder_ << std::endl;
    *out << "initialPhase_ = " << initialPhase_ << std::endl;
  }


  if (is_null(delta_)) {
    delta_ = createMember(stepper.get_x_space());
  }
  if (is_null(errWtVec_)) {
    errWtVec_ = createMember(stepper.get_x_space());
  }
  V_S(delta_.ptr(),zero); 

  setStepControlState_(BEFORE_FIRST_STEP);

  if (doTrace) {
    *out
      << "\nLeaving " << this->Teuchos::Describable::description()
      << "::initialize()...\n";
  }
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::getFirstTimeStep_(const StepperBase<Scalar>& stepper)
{
  
  TEST_FOR_EXCEPT(!(stepControlState_ == BEFORE_FIRST_STEP));

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar zero = ST::zero();

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool doTrace = (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH));
  Teuchos::OSTab ostab(out,1,"getFirstTimeStep_");

  const ImplicitBDFStepper<Scalar>& implicitBDFStepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
  const Thyra::VectorBase<Scalar>& xHistory0 = implicitBDFStepper.getxHistory(0);
  errWtVecCalc_->errWtVecSet(&*errWtVec_,xHistory0,relErrTol_,absErrTol_);

  // Choose initial step-size
  Scalar time_to_stop = stopTime_ - time_;
  Scalar currentTimeStep = ST::nan();
  if (stepSizeType_ == STEP_TYPE_FIXED) {
    currentTimeStep = hh_;
    //currentTimeStep = 0.1 * time_to_stop;
    //currentTimeStep = std::min(hh_, currentTimeStep);
  } else {
    // compute an initial step-size based on rate of change in the solution initially
    const Thyra::VectorBase<Scalar>& xHistory1 = implicitBDFStepper.getxHistory(1);
    Scalar ypnorm = wRMSNorm_(*errWtVec_,xHistory1);
    if (ypnorm > zero) { // time-dependent DAE
      currentTimeStep = std::min(h0_max_factor_*std::abs(time_to_stop),std::sqrt(2.0)/(h0_safety_*ypnorm));
    } else { // non-time-dependent DAE
      currentTimeStep = h0_max_factor_*std::abs(time_to_stop);
    }
    // choose std::min of user specified value and our value:
    if (hh_ > zero) {
      currentTimeStep = std::min(hh_, currentTimeStep);
    }
    // check for maximum step-size:
#ifdef RYTHMOS_DEBUG
      TEST_FOR_EXCEPT(ST::isnaninf(currentTimeStep));
#endif // RYTHMOS_DEBUG
    Scalar rh = std::abs(currentTimeStep)*h_max_inv_; 
    if (rh>1.0) {
      currentTimeStep = currentTimeStep/rh;
    }
  }
  hh_ = currentTimeStep;
  
  // Coefficient initialization 
  psi_[0] = hh_;
  cj_ = 1/psi_[0];

  if (doTrace) {
    *out << "\nhh_ = " << hh_ << std::endl;
    *out << "psi_[0] = " << psi_[0] << std::endl;
    *out << "cj_ = " << cj_ << std::endl;
  }
  
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setRequestedStepSize(
    const StepperBase<Scalar>& stepper
    ,const Scalar& stepSize
    ,const StepSizeType& stepSizeType
    )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT(!((stepControlState_ == UNINITIALIZED) || (stepControlState_ == BEFORE_FIRST_STEP) || (stepControlState_ == READY_FOR_NEXT_STEP) || (stepControlState_ == MID_STEP)));
  TEST_FOR_EXCEPTION(
      ((stepSizeType == STEP_TYPE_FIXED) && (stepSize == ST::zero())), 
      std::logic_error, 
      "Error, step size type == STEP_TYPE_FIXED, but requested step size == 0!\n"
      );
  bool didInitialization = false; 
  if (stepControlState_ == UNINITIALIZED) {
    initialize(stepper);
    didInitialization = true;
  }

  // errWtVecSet_ is called during initialize
  if (!didInitialization) {
    const ImplicitBDFStepper<Scalar>& implicitBDFStepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
    const Thyra::VectorBase<Scalar>& xHistory = implicitBDFStepper.getxHistory(0);
    errWtVecCalc_->errWtVecSet(&*errWtVec_,xHistory,relErrTol_,absErrTol_);
  }

  stepSizeType_ = stepSizeType;
  if (stepSizeType_ == STEP_TYPE_FIXED) {
    hh_ = stepSize;
    if (numberOfSteps_ == 0) {
      psi_[0] = hh_;
      cj_ = 1/psi_[0];
    }
  } else { // STEP_TYPE_VARIABLE
    if (stepSize != ST::zero()) {
      maxTimeStep_ = stepSize;
      h_max_inv_ = Scalar(ST::one()/maxTimeStep_);
    }
  }
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::nextStepSize(const StepperBase<Scalar>& stepper, Scalar* stepSize, StepSizeType* stepSizeType, int* order)
{
  TEST_FOR_EXCEPT(!((stepControlState_ == BEFORE_FIRST_STEP) || 
         (stepControlState_ == MID_STEP) ||  
         (stepControlState_ == READY_FOR_NEXT_STEP) )
        );
  if (stepControlState_ == BEFORE_FIRST_STEP) {
    getFirstTimeStep_(stepper);
  }
  if (stepControlState_ != MID_STEP) {
    this->updateCoeffs_();
  }
  if (stepSizeType_ == STEP_TYPE_VARIABLE) {
    if (hh_ > maxTimeStep_) {
      hh_ = maxTimeStep_;
    }
  }
  *stepSize = hh_;
  *stepSizeType = stepSizeType_;
  *order = currentOrder_;
  if (stepControlState_ != MID_STEP) {
    setStepControlState_(MID_STEP);
  }
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setCorrection(
     const StepperBase<Scalar>& stepper
    ,const RCP<const Thyra::VectorBase<Scalar> >& soln
    ,const RCP<const Thyra::VectorBase<Scalar> >& ee
    ,int solveStatus)
{
  TEST_FOR_EXCEPT(stepControlState_ != MID_STEP);
  TEST_FOR_EXCEPTION(is_null(ee), std::logic_error, "Error, ee == Teuchos::null!\n");
  ee_ = ee;
  newtonConvergenceStatus_ = solveStatus;
  setStepControlState_(AFTER_CORRECTION);
} 

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::completeStep(const StepperBase<Scalar>& stepper)
{
  TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  numberOfSteps_ ++;
  nef_ = 0;
  time_ += hh_;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"completeStep_");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "numberOfSteps_ = " << numberOfSteps_ << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
    *out << "time_ = " << time_ << std::endl;
  }
  
  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (stepSizeType_ == STEP_TYPE_VARIABLE);

  Scalar newTimeStep = hh_;
  Scalar rr = ST::one(); // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  int orderDiff = currentOrder_ - usedOrder_;
  usedOrder_ = currentOrder_;
  usedStep_ = hh_;
  if ((newOrder_ == currentOrder_-1) || (currentOrder_ == maxOrder_)) {
    // If we reduced our order or reached std::max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    initialPhase_ = false;
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "initialPhase_ = " << initialPhase_ << std::endl;
  }
  if (initialPhase_) {
    currentOrder_++;
    newTimeStep = h_phase0_incr_ * hh_;
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "currentOrder_ = " << currentOrder_ << std::endl;
      *out << "newTimeStep = " << newTimeStep << std::endl;
    }
  } else { // not in the initial phase of integration
    BDFactionFlag action = ACTION_UNSET;
    if (newOrder_ == currentOrder_-1) {
      action = ACTION_LOWER;
    } else if (newOrder_ == maxOrder_) {
      action = ACTION_MAINTAIN;
    } else if ((currentOrder_+1>=nscsco_) || (orderDiff == 1)) {
      // If we just raised the order last time then we won't raise it again
      // until we've taken currentOrder_+1 steps at order currentOrder_.
      action = ACTION_MAINTAIN;
    } else { // consider changing the order 
      const ImplicitBDFStepper<Scalar>& implicitBDFStepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
      const Thyra::VectorBase<Scalar>& xHistory = implicitBDFStepper.getxHistory(currentOrder_+1);
      V_StVpStV(delta_.ptr(),ST::one(),*ee_,Scalar(-ST::one()),xHistory);
      Tkp1_ = wRMSNorm_(*errWtVec_,*delta_);
      Ekp1_ = Tkp1_/(currentOrder_+2);
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "delta_ = " << std::endl;
        delta_->describe(*out,verbLevel);
        *out << "Tkp1_ = ||delta_||_WRMS = " << Tkp1_ << std::endl;
        *out << "Ekp1_ = " << Ekp1_ << std::endl;
      }
      if (currentOrder_ == 1) {
        if (Tkp1_ >= Tkp1_Tk_safety_ * Tk_) {
          action = ACTION_MAINTAIN;
        } else {
          action = ACTION_RAISE;
        }
      } else {
        if (Tkm1_ <= std::min(Tk_,Tkp1_)) {
          action = ACTION_LOWER;
        } else if (Tkp1_ >= Tk_) {
          action = ACTION_MAINTAIN;
        } else {
          action = ACTION_RAISE;
        }
      }
    }
    if (currentOrder_ < minOrder_) {
      action = ACTION_RAISE;
    } else if ( (currentOrder_ == minOrder_) && (action == ACTION_LOWER) ) {
      action = ACTION_MAINTAIN;
    }
    if (action == ACTION_RAISE) {
      currentOrder_++;
      Est_ = Ekp1_;
    } else if (action == ACTION_LOWER) {
      currentOrder_--;
      Est_ = Ekm1_;
    }
    newTimeStep = hh_;
    rr = pow(r_safety_*Est_+r_fudge_,-1.0/(currentOrder_+1.0));
    if (rr >= r_hincr_test_) {
      rr = r_hincr_;
      newTimeStep = rr*hh_;
    } else if (rr <= 1) {
      rr = std::max(r_min_,std::min(r_max_,rr));
      newTimeStep = rr*hh_;
    }
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Est_ = " << Est_ << std::endl;
      *out << "rr  = " << rr << std::endl;
      *out << "newTimeStep = " << newTimeStep << std::endl;
    }
  }
  
  if (time_ < stopTime_) {
    // If the step needs to be adjusted:
    if (adjustStep) {
      newTimeStep = std::max(newTimeStep, minTimeStep_);
      newTimeStep = std::min(newTimeStep, maxTimeStep_);

      Scalar nextTimePt = time_ + newTimeStep;

      if (nextTimePt > stopTime_) {
        nextTimePt  = stopTime_;
        newTimeStep = stopTime_ - time_;
      }

      hh_ = newTimeStep;

    } else { // if time step is constant for this step:
      Scalar nextTimePt = time_ + hh_;

      if (nextTimePt > stopTime_) {
        nextTimePt      = stopTime_;
        hh_ = stopTime_ - time_;
      }
    }
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "hh_ = " << hh_ << std::endl;
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
  setStepControlState_(READY_FOR_NEXT_STEP);
}

template<class Scalar>
AttemptedStepStatusFlag ImplicitBDFStepperStepControl<Scalar>::rejectStep(const StepperBase<Scalar>& stepper)
{
  TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);

  using Teuchos::as;

  // This routine puts its output in newTimeStep and newOrder

  // This routine changes the following variables:
  //    initialPhase, nef, psi, newTimeStep,
  //    newOrder, currentOrder_, currentTimeStep, dsDae.xHistory,
  //    dsDae.qHistory, nextTimePt, 
  //    currentTimeStepSum, nextTimePt

  // This routine reads but does not change the following variables:
  //    r_factor, r_safety, Est_, r_fudge_, r_min_, r_max_,
  //    minTimeStep_, maxTimeStep_, time, stopTime_ 

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (stepSizeType_ == STEP_TYPE_VARIABLE);

  Scalar newTimeStep = hh_;
  Scalar rr = 1.0; // step size ratio = new step / old step
  nef_++;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"rejectStep_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "adjustStep = " << adjustStep << std::endl;
    *out << "nef_ = " << nef_ << std::endl;
  }
  if (nef_ >= max_LET_fail_)  {
    TEST_FOR_EXCEPTION(nef_ >= max_LET_fail_, std::logic_error, "Error, maximum number of local error test failures.\n");
  }
  initialPhase_ = false;
  if (adjustStep) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "initialPhase_ = " << initialPhase_ << std::endl;
    }
    for (int i=1;i<=currentOrder_;++i) {
      psi_[i-1] = psi_[i] - hh_;
    }

    if ((newtonConvergenceStatus_ < 0)) {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = r_min_;
      newTimeStep = rr * hh_;

      if (nef_ > 2) {
        newOrder_ = 1;//consistent with block below.
      }
    } else {
      // 03/11/04 tscoffe:  Here is the block for choosing order & 
      // step-size when the local error test FAILS (but Newton 
      // succeeded). 
      if (nef_ == 1) { // first local error test failure
        rr = r_factor_*pow(r_safety_*Est_+r_fudge_,-1.0/(newOrder_+1.0));
        rr = std::max(r_min_,std::min(r_max_,rr));
        newTimeStep = rr * hh_;
      } else if (nef_ == 2) { // second failure
        rr = r_min_;
        newTimeStep = rr * hh_;
      } else if (nef_ > 2) { // third and later failures
        newOrder_ = 1;
        rr = r_min_;
        newTimeStep = rr * hh_;
      }
    }
    if (newOrder_ >= minOrder_) {
      currentOrder_ = newOrder_;
    }
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "rr = " << rr << std::endl;
      *out << "newOrder_ = " << newOrder_ << std::endl;
      *out << "currentOrder_ = " << currentOrder_ << std::endl;
    }
    if (numberOfSteps_ == 0) { // still first step
      psi_[0] = newTimeStep;
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "numberOfSteps_ == 0:" << std::endl;
        *out << "psi_[0] = " << psi_[0] << std::endl;
      }
    }
  } else if (!adjustStep) {
    if ( as<int>(verbLevel) != as<int>(Teuchos::VERB_NONE) ) {
      *out << "Rythmos_ImplicitBDFStepperStepControl::rejectStep(...):  "
          << "Warning:  Local error test failed with constant step-size."
          << std::endl;
    }
  }

  AttemptedStepStatusFlag return_status = PREDICT_AGAIN;

  // If the step needs to be adjusted:
  if (adjustStep) {
    newTimeStep = std::max(newTimeStep, minTimeStep_);
    newTimeStep = std::min(newTimeStep, maxTimeStep_);

    Scalar nextTimePt = time_ + newTimeStep;

    if (nextTimePt > stopTime_) {
      nextTimePt  = stopTime_;
      newTimeStep = stopTime_ - time_;
    }

    hh_ = newTimeStep;
  
  } else { // if time step is constant for this step:
    Scalar nextTimePt = time_ + hh_;

    if (nextTimePt > stopTime_) {
      nextTimePt      = stopTime_;
      hh_ = stopTime_ - time_;
    }
    return_status = CONTINUE_ANYWAY;
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "hh_ = " << hh_ << std::endl;
  }

  if (return_status == PREDICT_AGAIN) {
    setStepControlState_(READY_FOR_NEXT_STEP);
  } else if (return_status == CONTINUE_ANYWAY) {
    // do nothing, as we'll call completeStep which must be in AFTER_CORRECTION state.
  }

  return(return_status);
}

template<class Scalar>
Scalar ImplicitBDFStepperStepControl<Scalar>::checkReduceOrder_(const StepperBase<Scalar>& stepper)
{
  TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);
  TEST_FOR_EXCEPT(checkReduceOrderCalled_ == true);

  using Teuchos::as;

  const ImplicitBDFStepper<Scalar>& implicitBDFStepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);

  // This routine puts its output in newOrder_
  
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"checkReduceOrder_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "sigma_[" << currentOrder_ << "] = " << sigma_[currentOrder_] << std::endl;
    *out << "ee_ = " << std::endl;
    ee_->describe(*out,verbLevel);
    *out << "errWtVec_ = " << std::endl;
    errWtVec_->describe(*out,verbLevel);
  }

  Scalar enorm = wRMSNorm_(*errWtVec_,*ee_);
  Ek_ = sigma_[currentOrder_]*enorm;
  Tk_ = Scalar(currentOrder_+1)*Ek_;
  Est_ = Ek_;
  newOrder_ = currentOrder_;
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
    *out << "Ek_ = " << Ek_ << std::endl;
    *out << "Tk_ = " << Tk_ << std::endl;
    *out << "enorm = " << enorm << std::endl;
  }
  if (currentOrder_>1) {
    const Thyra::VectorBase<Scalar>& xHistoryCur = implicitBDFStepper.getxHistory(currentOrder_);
    V_VpV(delta_.ptr(),xHistoryCur,*ee_);
    Ekm1_ = sigma_[currentOrder_-1]*wRMSNorm_(*errWtVec_,*delta_);
    Tkm1_ = currentOrder_*Ekm1_;
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Ekm1_ = " << Ekm1_ << std::endl;
      *out << "Tkm1_ = " << Tkm1_ << std::endl;
    }
    if (currentOrder_>2) {
      const Thyra::VectorBase<Scalar>& xHistoryPrev = implicitBDFStepper.getxHistory(currentOrder_-1);
      Vp_V(delta_.ptr(),xHistoryPrev);
      Ekm2_ = sigma_[currentOrder_-2]*wRMSNorm_(*errWtVec_,*delta_);
      Tkm2_ = (currentOrder_-1)*Ekm2_;
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Ekm2_ = " << Ekm2_ << std::endl;
        *out << "Tkm2_ = " << Tkm2_ << std::endl;
      }
      if (std::max(Tkm1_,Tkm2_)<=Tk_) {
        newOrder_--;
        Est_ = Ekm1_;
      }
    }
    else if (Tkm1_ <= Tkm1_Tk_safety_ * Tk_) {
      newOrder_--;
      Est_ = Ekm1_;
    }
  }
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Est_ = " << Est_ << std::endl;
    *out << "newOrder_= " << newOrder_ << std::endl;
  }
  checkReduceOrderCalled_ = true;
  return(enorm);
}

template<class Scalar>
bool ImplicitBDFStepperStepControl<Scalar>::acceptStep(const StepperBase<Scalar>& stepper, Scalar* LETValue)
{
  TEST_FOR_EXCEPT(stepControlState_ != AFTER_CORRECTION);
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"acceptStep");

  bool return_status = false;
  Scalar enorm = checkReduceOrder_(stepper);
  Scalar LET = ck_*enorm;

  if (failStepIfNonlinearSolveFails_ && (newtonConvergenceStatus_ < 0) )
    return false;

  if (LETValue) {
    *LETValue = LET;
  }
  if (LET < ST::one()) {
    return_status = true;
  }
  if ( Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "ck_ = " << ck_ << std::endl;
    *out << "enorm = " << enorm << std::endl;
    *out << "Local Truncation Error Check: (ck*enorm) < 1:  (" << LET << ") <?= 1" << std::endl;
  }
  return(return_status);
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::describe(
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
    out << "time_ = " << time_ << std::endl;
    out << "hh_ = " << hh_ << std::endl;
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
void ImplicitBDFStepperStepControl<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEST_FOR_EXCEPT(paramList == Teuchos::null);
  paramList->validateParameters(*this->getValidParameters(),0);
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);

  minOrder_ = parameterList_->get("minOrder",int(1)); // minimum order
  TEST_FOR_EXCEPTION(
      !((1 <= minOrder_) && (minOrder_ <= 5)), std::logic_error,
      "Error, minOrder_ = " << minOrder_ << " is not in range [1,5]!\n"
      );
  maxOrder_ = parameterList_->get("maxOrder",int(5)); // maximum order
  TEST_FOR_EXCEPTION(
      !((1 <= maxOrder_) && (maxOrder_ <= 5)), std::logic_error,
      "Error, maxOrder_ = " << maxOrder_ << " is not in range [1,5]!\n"
      );

  relErrTol_ = parameterList_->get( "relErrTol", Scalar(1.0e-4) );
  absErrTol_ = parameterList_->get( "absErrTol", Scalar(1.0e-6) );
  bool constantStepSize = parameterList_->get( "constantStepSize", false );
  stopTime_ = parameterList_->get( "stopTime", Scalar(1.0) );
  
  if (constantStepSize == true) {
    stepSizeType_ = STEP_TYPE_FIXED;
  } else {
    stepSizeType_ = STEP_TYPE_VARIABLE;
  }

  failStepIfNonlinearSolveFails_ = 
    parameterList_->get( "failStepIfNonlinearSolveFails", false );

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setParameterList");
  out->precision(15);

  setDefaultMagicNumbers_(parameterList_->sublist("magicNumbers"));

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "minOrder_ = " << minOrder_ << std::endl;
    *out << "maxOrder_ = " << maxOrder_ << std::endl;
    *out << "relErrTol  = " << relErrTol_  << std::endl;
    *out << "absErrTol  = " << absErrTol_  << std::endl;
    *out << "stepSizeType = " << stepSizeType_  << std::endl;
    *out << "stopTime_  = " << stopTime_  << std::endl;
    *out << "failStepIfNonlinearSolveFails_ = " 
	 << failStepIfNonlinearSolveFails_  << std::endl;
  }

}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setDefaultMagicNumbers_(
  Teuchos::ParameterList &magicNumberList)
{

  using Teuchos::as;

  // Magic Number Defaults:
  h0_safety_      = magicNumberList.get( "h0_safety",      Scalar(2.0)     );
  h0_max_factor_  = magicNumberList.get( "h0_max_factor",  Scalar(0.001)   );
  h_phase0_incr_  = magicNumberList.get( "h_phase0_incr",  Scalar(2.0)     );
  h_max_inv_      = magicNumberList.get( "h_max_inv",      Scalar(0.0)     );
  Tkm1_Tk_safety_ = magicNumberList.get( "Tkm1_Tk_safety", Scalar(2.0)     );
  Tkp1_Tk_safety_ = magicNumberList.get( "Tkp1_Tk_safety", Scalar(0.5)     );
  r_factor_       = magicNumberList.get( "r_factor",       Scalar(0.9)     );
  r_safety_       = magicNumberList.get( "r_safety",       Scalar(2.0)     );
  r_fudge_        = magicNumberList.get( "r_fudge",        Scalar(0.0001)  );
  r_min_          = magicNumberList.get( "r_min",          Scalar(0.125)   );
  r_max_          = magicNumberList.get( "r_max",          Scalar(0.9)     );
  r_hincr_test_   = magicNumberList.get( "r_hincr_test",   Scalar(2.0)     );
  r_hincr_        = magicNumberList.get( "r_hincr",        Scalar(2.0)     );
  max_LET_fail_   = magicNumberList.get( "max_LET_fail",   int(15)         );
  minTimeStep_    = magicNumberList.get( "minTimeStep",    Scalar(0.0)     );
  maxTimeStep_    = magicNumberList.get( "maxTimeStep",    Scalar(10.0)    ); 

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setDefaultMagicNumbers_");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "h0_safety_ = " << h0_safety_ << std::endl;
    *out << "h0_max_factor_ = " << h0_max_factor_ << std::endl;
    *out << "h_phase0_incr_ = " << h_phase0_incr_ << std::endl;
    *out << "h_max_inv_ = " << h_max_inv_ << std::endl;
    *out << "Tkm1_Tk_safety_ = " << Tkm1_Tk_safety_  << std::endl;
    *out << "Tkp1_Tk_safety_ = " << Tkp1_Tk_safety_ << std::endl;
    *out << "r_factor_ = " << r_factor_ << std::endl;
    *out << "r_safety_ = " << r_safety_ << std::endl;
    *out << "r_fudge_ = " << r_fudge_ << std::endl;
    *out << "r_min_ = " << r_min_ << std::endl;
    *out << "r_max_ = " << r_max_ << std::endl;
    *out << "r_hincr_test_ = " << r_hincr_test_ << std::endl;
    *out << "r_hincr_ = " << r_hincr_ << std::endl;
    *out << "max_LET_fail_ = " << max_LET_fail_ << std::endl;
    *out << "minTimeStep_ = " << minTimeStep_ << std::endl;
    *out << "maxTimeStep_ = " << maxTimeStep_ << std::endl;
  }

}

template<class Scalar>
RCP<const Teuchos::ParameterList>
ImplicitBDFStepperStepControl<Scalar>::getValidParameters() const
{

  static RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();

    pl->set<int>   ( "minOrder",         1,
        "lower limit of order selection, guaranteed"
        );
    pl->set<int>   ( "maxOrder",         5,
        "upper limit of order selection, does not guarantee this order"        
        );
    pl->set<Scalar>( "relErrTol",        Scalar(1.0e-4) );
    pl->set<Scalar>( "absErrTol",        Scalar(1.0e-6) );
    pl->set<bool>  ( "constantStepSize", false          );
    pl->set<Scalar>( "stopTime",         Scalar(10.0)   );
    pl->set<bool>("failStepIfNonlinearSolveFails", false,
		  "Power user command. Will force the function acceptStep() to return false ieven if the LET is acceptable.  Used to run with loose tolerances but enforce a correct nonlinear solution to the step.");

    Teuchos::ParameterList
      &magicNumberList = pl->sublist("magicNumbers", 
          false,
          "These are knobs in the algorithm that have been set to reasonable values using lots of testing and heuristics and some theory."
          );
    magicNumberList.set<Scalar>( "h0_safety",      Scalar(2.0)     );
    magicNumberList.set<Scalar>( "h0_max_factor",  Scalar(0.001)   );
    magicNumberList.set<Scalar>( "h_phase0_incr",  Scalar(2.0),
        "initial ramp-up in variable mode (stepSize multiplier) "     
        );
    magicNumberList.set<Scalar>( "h_max_inv",      Scalar(0.0)     );
    magicNumberList.set<Scalar>( "Tkm1_Tk_safety", Scalar(2.0)     );
    magicNumberList.set<Scalar>( "Tkp1_Tk_safety", Scalar(0.5)     );
    magicNumberList.set<Scalar>( "r_factor",       Scalar(0.9),
        "used in rejectStep:  time step ratio multiplier"
        );
    magicNumberList.set<Scalar>( "r_safety",       Scalar(2.0),
        "local error multiplier as part of time step ratio calculation"
        );
    magicNumberList.set<Scalar>( "r_fudge",        Scalar(0.0001),
        "local error addition as part of time step ratio calculation"
        );
    magicNumberList.set<Scalar>( "r_min",          Scalar(0.125),
        "used in rejectStep:  how much to cut step and lower bound for time step ratio"   
        );
    magicNumberList.set<Scalar>( "r_max",          Scalar(0.9),
        "upper bound for time step ratio"
        );
    magicNumberList.set<Scalar>( "r_hincr_test",   Scalar(2.0),     
        "used in completeStep:  if time step ratio > this then set time step ratio to r_hincr"
        );
    magicNumberList.set<Scalar>( "r_hincr",        Scalar(2.0),
        "used in completeStep:  limit on time step ratio increases, not checked by r_max"
        );
    magicNumberList.set<int>   ( "max_LET_fail",   int(15),
        "Max number of rejected steps"
        );
    magicNumberList.set<Scalar>( "minTimeStep",    Scalar(0.0),
        "bound on smallest time step in variable mode."     
        );
    magicNumberList.set<Scalar>( "maxTimeStep",    Scalar(10.0),
        "bound on largest time step in variable mode."    
        ); 

    Teuchos::setupVerboseObjectSublist(&*pl);

    validPL = pl;

  }

  return (validPL);
  
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperStepControl<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
ImplicitBDFStepperStepControl<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setStepControlData(const StepperBase<Scalar>& stepper)
{
  if (stepControlState_ == UNINITIALIZED) {
    initialize(stepper);
  }
  const ImplicitBDFStepper<Scalar>& bdfstepper = Teuchos::dyn_cast<const ImplicitBDFStepper<Scalar> >(stepper);
  int desiredOrder = bdfstepper.getOrder();
  TEST_FOR_EXCEPT(!((1 <= desiredOrder) && (desiredOrder <= maxOrder_)));
  if (stepControlState_ == BEFORE_FIRST_STEP) {
    TEST_FOR_EXCEPTION(
        desiredOrder > 1, 
        std::logic_error, 
        "Error, this ImplicitBDF stepper has not taken a step yet, so it cannot take a step of order " << desiredOrder << " > 1!\n"
        );
  }
  TEST_FOR_EXCEPT(!(desiredOrder <= currentOrder_+1));
  currentOrder_ = desiredOrder;

  using Teuchos::as;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::OSTab ostab(out,1,"setStepControlData");

  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) ) {
    *out << "currentOrder_ = " << currentOrder_ << std::endl;
  }
}

template<class Scalar>
bool ImplicitBDFStepperStepControl<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
ImplicitBDFStepperStepControl<Scalar>::cloneStepControlStrategyAlgorithm() const
{

  RCP<ImplicitBDFStepperStepControl<Scalar> > stepControl = rcp(new ImplicitBDFStepperStepControl<Scalar>());

  if (!is_null(parameterList_)) {
    stepControl->setParameterList(parameterList_);
  }

  return stepControl;
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::setErrWtVecCalc(const RCP<ErrWtVecCalcBase<Scalar> >& errWtVecCalc)
{
  TEST_FOR_EXCEPT(is_null(errWtVecCalc));
  errWtVecCalc_ = errWtVecCalc;
}

template<class Scalar>
RCP<const ErrWtVecCalcBase<Scalar> > ImplicitBDFStepperStepControl<Scalar>::getErrWtVecCalc() const
{
  return(errWtVecCalc_);
}

template<class Scalar>
Scalar ImplicitBDFStepperStepControl<Scalar>::wRMSNorm_(
    const Thyra::VectorBase<Scalar>& weight, 
    const Thyra::VectorBase<Scalar>& vector) const
{
  return(norm_2(weight,vector));
}

template<class Scalar>
void ImplicitBDFStepperStepControl<Scalar>::defaultInitializeAllData_()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar zero = ST::zero();
  Scalar mone = Scalar(-ST::one());

  stepControlState_ = UNINITIALIZED;
  hh_ = zero;
  numberOfSteps_ = 0;
  stepSizeType_ = STEP_TYPE_VARIABLE;

  minOrder_ = -1;
  maxOrder_ = -1;
  nef_ = 0;
  midStep_ = false;
  checkReduceOrderCalled_ = false;
  time_ = -std::numeric_limits<Scalar>::max();
  relErrTol_ = mone;
  absErrTol_ = mone;
  usedStep_ = mone;
  currentOrder_ = 1;
  usedOrder_ = -1;
  nscsco_ = -1;
  alpha_s_ = mone;
  alpha_0_ = mone;
  cj_ = mone;
  ck_ = mone;
  ck_enorm_ = mone;
  constantStepSize_ = false;
  Ek_ = mone;
  Ekm1_ = mone;
  Ekm2_ = mone;
  Ekp1_ = mone;
  Est_ = mone;
  Tk_ = mone;
  Tkm1_ = mone;
  Tkm2_ = mone;
  Tkp1_ = mone;
  newOrder_ = -1;
  oldOrder_ = -1;
  initialPhase_ = false;
  stopTime_ = mone;
  h0_safety_ = mone;
  h0_max_factor_ = mone;
  h_phase0_incr_ = mone;
  h_max_inv_ = mone;
  Tkm1_Tk_safety_ = mone;
  Tkp1_Tk_safety_ = mone;
  r_factor_ = mone;
  r_safety_ = mone;
  r_fudge_ = mone;
  r_min_ = mone;
  r_max_ = mone;
  r_hincr_test_ = mone;
  r_hincr_ = mone;
  max_LET_fail_ = -1;
  minTimeStep_ = mone;
  maxTimeStep_ = mone;
  newtonConvergenceStatus_ = -1;
}

template<class Scalar>
int ImplicitBDFStepperStepControl<Scalar>::getMinOrder() const
{
  TEST_FOR_EXCEPTION(
      stepControlState_ == UNINITIALIZED, std::logic_error,
      "Error, attempting to call getMinOrder before intiialization!\n"
      );
  return(minOrder_);
}

template<class Scalar>
int ImplicitBDFStepperStepControl<Scalar>::getMaxOrder() const
{
  TEST_FOR_EXCEPTION(
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

#define RYTHMOS_IMPLICITBDF_STEPPER_STEPCONTROL_INSTANT(SCALAR) \
  template class ImplicitBDFStepperStepControl< SCALAR >; 


} // namespace Rythmos
#endif // Rythmos_IMPLICITBDF_STEPPER_STEP_CONTROL_DEF_H

