// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SolutionState_impl_hpp
#define Tempus_SolutionState_impl_hpp

#include "Thyra_VectorStdOps.hpp"

namespace Tempus {


template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
  : x_              (x),
    x_nc_           (x),
    xdot_           (xdot),
    xdot_nc_        (xdot),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (xdotdot),
    stepperState_   (stepperState),
    stepperState_nc_(stepperState),
    physicsState_   (physicsState),
    physicsState_nc_(physicsState)
{
  metaData_nc_ = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
  metaData_    = metaData_nc_;
  if (stepperState_nc_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_nc_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
    physicsState_    = physicsState_nc_;
  }
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<const StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<const PhysicsState<Scalar> >& physicsState)
  : x_              (x),
    x_nc_           (Teuchos::null),
    xdot_           (xdot),
    xdot_nc_        (Teuchos::null),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (Teuchos::null),
    stepperState_   (Teuchos::null),
    stepperState_nc_(Teuchos::null),
    physicsState_   (Teuchos::null),
    physicsState_nc_(Teuchos::null)
{
  metaData_nc_ = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
  metaData_    = metaData_nc_;
  stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
  stepperState_    = stepperState_nc_;
  physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
  physicsState_    = physicsState_nc_;
}


template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
  : metaData_       (metaData),
    metaData_nc_    (metaData),
    x_              (x),
    x_nc_           (x),
    xdot_           (xdot),
    xdot_nc_        (xdot),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (xdotdot),
    stepperState_   (stepperState),
    stepperState_nc_(stepperState),
    physicsState_   (physicsState),
    physicsState_nc_(physicsState)
{
  if (stepperState_nc_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_nc_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
    physicsState_    = physicsState_nc_;
  }
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<const SolutionStateMetaData<Scalar> > metaData,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<const StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<const PhysicsState<Scalar> >& physicsState)
  : metaData_       (metaData),
    metaData_nc_    (Teuchos::null),
    x_              (x),
    x_nc_           (Teuchos::null),
    xdot_           (xdot),
    xdot_nc_        (Teuchos::null),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (Teuchos::null),
    stepperState_   (stepperState),
    stepperState_nc_(Teuchos::null),
    physicsState_   (physicsState),
    physicsState_nc_(Teuchos::null)
{
  if (stepperState_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
    physicsState_    = physicsState_nc_;
  }
}


template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Scalar time,
  const Scalar dt,
  const int    iStep,
  const Scalar errorAbs,
  const Scalar errorRel,
  const int    order,
  const int    nFailures,
  const int    nRunningFailures,
  const int    nConsecutiveFailures,
  const Status solutionStatus,
  const bool   output,
  const bool   outputScreen,
  const bool   isSynced,
  const bool   isInterpolated,
  const Scalar accuracy,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
  : x_              (x),
    x_nc_           (x),
    xdot_           (xdot),
    xdot_nc_        (xdot),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (xdotdot),
    stepperState_   (stepperState),
    stepperState_nc_(stepperState),
    physicsState_   (physicsState),
    physicsState_nc_(physicsState)
{
  metaData_nc_ =
    Teuchos::rcp(new SolutionStateMetaData<Scalar> (time,
                                                    iStep,
                                                    dt,
                                                    errorAbs,
                                                    errorRel,
                                                    order,
                                                    nFailures,
                                                    nRunningFailures,
                                                    nConsecutiveFailures,
                                                    solutionStatus,
                                                    output,
                                                    outputScreen,
                                                    isSynced,
                                                    isInterpolated,
                                                    accuracy));
  metaData_ = metaData_nc_;

  if (stepperState_nc_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_nc_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
    physicsState_ = physicsState_nc_;
  }
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Scalar time,
  const Scalar dt,
  const int    iStep,
  const Scalar errorAbs,
  const Scalar errorRel,
  const int    order,
  const int    nFailures,
  const int    nRunningFailures,
  const int    nConsecutiveFailures,
  const Status solutionStatus,
  const bool   output,
  const bool   outputScreen,
  const bool   isSynced,
  const bool   isInterpolated,
  const Scalar accuracy,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<const StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<const PhysicsState<Scalar> >& physicsState)
  : x_              (x),
    x_nc_           (Teuchos::null),
    xdot_           (xdot),
    xdot_nc_        (Teuchos::null),
    xdotdot_        (xdotdot),
    xdotdot_nc_     (Teuchos::null),
    stepperState_   (stepperState),
    stepperState_nc_(Teuchos::null),
    physicsState_   (physicsState),
    physicsState_nc_(Teuchos::null)
{
  metaData_nc_ =
    Teuchos::rcp(new SolutionStateMetaData<Scalar> (time,
                                                    iStep,
                                                    dt,
                                                    errorAbs,
                                                    errorRel,
                                                    order,
                                                    nFailures,
                                                    nRunningFailures,
                                                    nConsecutiveFailures,
                                                    solutionStatus,
                                                    output,
                                                    outputScreen,
                                                    isSynced,
                                                    isInterpolated,
                                                    accuracy));
  metaData_ = metaData_nc_;

  if (stepperState_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ());
    physicsState_ = physicsState_nc_;
  }
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  const Teuchos::RCP<StepperState<Scalar> >& stepperState,
  const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::rcp_const_cast;

  metaData_nc_ = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
  metaData_nc_->setSolutionStatus(Status::PASSED);
  metaData_ = metaData_nc_;

  MEB::InArgs<Scalar> inArgs = model->getNominalValues();

  // The solution vector, x, is required (usually).
  x_nc_ = rcp_const_cast<Thyra::VectorBase<Scalar> > (inArgs.get_x());
  x_    = x_nc_;

  // The solution derivative, xdot, can be optional provided, based on
  // application needs.  Here we will base it on "supports" IN_ARG_x_dot.
  // Depending on the stepper used, a temporary xdot vector may be created
  // within the Stepper, but not moved to the SolutionState.
  if (inArgs.supports(MEB::IN_ARG_x_dot)) {
    xdot_nc_ = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
    xdot_    = xdot_nc_;
  } else {
    xdot_nc_ = Teuchos::null;
    xdot_    = xdot_nc_;
  }

  // Similar as xdot.
  if (inArgs.supports(MEB::IN_ARG_x_dot_dot)) {
    xdotdot_nc_ =
      rcp_const_cast<Thyra::VectorBase<Scalar> > (inArgs.get_x_dot_dot());
    xdotdot_    = xdotdot_nc_;
  } else {
    xdotdot_nc_ = Teuchos::null;
    xdotdot_    = xdotdot_nc_;
  }

  if (stepperState_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar> ()); // Use default
    stepperState_    = stepperState_nc_;
  } else {
    stepperState_nc_ = stepperState;
    stepperState_    = stepperState;
  }

  if (physicsState_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar> ()); // Use default
    physicsState_    = physicsState_nc_;
  } else {
    physicsState_nc_ = physicsState;
    physicsState_    = physicsState;
  }
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(const SolutionState<Scalar>& ss_)
  :metaData_       (ss_.metaData_),
   metaData_nc_    (ss_.metaData_nc_),
   x_              (ss_.x_),
   x_nc_           (ss_.x_nc_),
   xdot_           (ss_.xdot_),
   xdot_nc_        (ss_.xdot_nc_),
   xdotdot_        (ss_.xdotdot_),
   xdotdot_nc_     (ss_.xdotdot_nc_),
   stepperState_   (ss_.stepperState_),
   stepperState_nc_(ss_.stepperState_nc_),
   physicsState_   (ss_.physicsState_),
   physicsState_nc_(ss_.physicsState_nc_)
{}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> > SolutionState<Scalar>::clone() const
{
  using Teuchos::RCP;

  RCP<SolutionStateMetaData<Scalar> > metaData_out;
  if (!Teuchos::is_null(metaData_)) metaData_out = metaData_->clone();

  RCP<Thyra::VectorBase<Scalar> > x_out;
  if (!Teuchos::is_null(x_)) x_out = x_->clone_v();

  RCP<Thyra::VectorBase<Scalar> > xdot_out;
  if (!Teuchos::is_null(xdot_)) xdot_out = xdot_->clone_v();

  RCP<Thyra::VectorBase<Scalar> > xdotdot_out;
  if (!Teuchos::is_null(xdotdot_)) xdotdot_out = xdotdot_->clone_v();

  RCP<StepperState<Scalar> > sS_out;
  if (!Teuchos::is_null(stepperState_)) sS_out=stepperState_->clone();

  RCP<PhysicsState<Scalar> > pS_out;
  if (!Teuchos::is_null(physicsState_)) pS_out=physicsState_->clone();

  RCP<SolutionState<Scalar> > ss_out = Teuchos::rcp(new SolutionState<Scalar> (
    metaData_out, x_out, xdot_out, xdotdot_out, sS_out, pS_out));

  return ss_out;
}


template<class Scalar>
void SolutionState<Scalar>::
copy(const Teuchos::RCP<const SolutionState<Scalar> >& ss)
{
  metaData_nc_->copy(ss->metaData_);
  this->copySolutionStepperState(ss);
}


template<class Scalar>
void SolutionState<Scalar>::
copySolutionStepperState(const Teuchos::RCP<const SolutionState<Scalar> >& ss)
{
  Thyra::V_V(x_nc_.ptr(),       *(ss->x_));
  if (ss->xdot_ == Teuchos::null) xdot_nc_ = Teuchos::null;
  else Thyra::V_V(xdot_nc_.ptr(),    *(ss->xdot_));
  if (ss->xdotdot_ == Teuchos::null) xdotdot_nc_ = Teuchos::null;
  else Thyra::V_V(xdotdot_nc_.ptr(), *(ss->xdotdot_));
  stepperState_nc_->copy(ss->stepperState_);
  physicsState_nc_->copy(ss->physicsState_);
}


template<class Scalar>
bool SolutionState<Scalar>::operator< (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() < ss.metaData_->getTime());
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() <= ss.metaData_->getTime());
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const Scalar& t) const
{
  return (this->metaData_->getTime() < t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const Scalar& t) const
{
  return (this->metaData_->getTime() <= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() > ss.metaData_->getTime());
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() >= ss.metaData_->getTime());
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const Scalar& t) const
{
  return (this->metaData_->getTime() > t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const Scalar& t) const
{
  return (this->metaData_->getTime() >= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() == ss.metaData_->getTime());
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const Scalar& t) const
{
  return (this->metaData_->getTime() == t);
}

template<class Scalar>
std::string SolutionState<Scalar>::description() const
{
  std::string name = "Tempus::SolutionState";
  return (name);
}

template<class Scalar>
void SolutionState<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  if (verbLevel == Teuchos::VERB_EXTREME) {
    out << description() << "::describe:" << std::endl
        << "metaData = " << std::endl;
        metaData_->describe(out,verbLevel);
    out << "x = " << std::endl;
    x_->describe(out,verbLevel);
    if (xdot_ != Teuchos::null) {
      out << "xdot_ = " << std::endl;
      xdot_->describe(out,verbLevel);
    }
    if (xdotdot_ != Teuchos::null) {
      out << "xdotdot = " << std::endl;
      xdotdot_->describe(out,verbLevel);
    }
    if (stepperState_ != Teuchos::null) {
      out << "stepperState = " << std::endl;
      stepperState_->describe(out,verbLevel);
    }
    if (physicsState_ != Teuchos::null) {
      out << "stepperState = " << std::endl;
      physicsState_->describe(out,verbLevel);
    }
  }
}

} // namespace Tempus
#endif // Tempus_SolutionState_impl_hpp
