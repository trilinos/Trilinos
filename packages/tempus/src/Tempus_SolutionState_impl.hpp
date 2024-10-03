//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_SolutionState_impl_hpp
#define Tempus_SolutionState_impl_hpp

#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
SolutionState<Scalar>::SolutionState()
  : x_(Teuchos::null),
    x_nc_(Teuchos::null),
    xdot_(Teuchos::null),
    xdot_nc_(Teuchos::null),
    xdotdot_(Teuchos::null),
    xdotdot_nc_(Teuchos::null),
    stepperState_(Teuchos::null),
    stepperState_nc_(Teuchos::null),
    physicsState_(Teuchos::null),
    physicsState_nc_(Teuchos::null)
{
  metaData_nc_     = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
  metaData_        = metaData_nc_;
  stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
  stepperState_    = stepperState_nc_;
  physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar>());
  physicsState_    = physicsState_nc_;
}

template <class Scalar>
SolutionState<Scalar>::SolutionState(
    const Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<StepperState<Scalar> >& stepperState,
    const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
  : metaData_(metaData),
    metaData_nc_(metaData),
    x_(x),
    x_nc_(x),
    xdot_(xdot),
    xdot_nc_(xdot),
    xdotdot_(xdotdot),
    xdotdot_nc_(xdotdot),
    stepperState_(stepperState),
    stepperState_nc_(stepperState),
    physicsState_(physicsState),
    physicsState_nc_(physicsState)
{
  if (stepperState_nc_ == Teuchos::null) {
    stepperState_nc_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
    stepperState_    = stepperState_nc_;
  }
  if (physicsState_nc_ == Teuchos::null) {
    physicsState_nc_ = Teuchos::rcp(new PhysicsState<Scalar>());
    physicsState_    = physicsState_nc_;
  }
}

template <class Scalar>
SolutionState<Scalar>::SolutionState(
    const Teuchos::RCP<const SolutionStateMetaData<Scalar> > metaData,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot,
    const Teuchos::RCP<const StepperState<Scalar> >& stepperState,
    const Teuchos::RCP<const PhysicsState<Scalar> >& physicsState)
  : metaData_(metaData),
    metaData_nc_(Teuchos::null),
    x_(x),
    x_nc_(Teuchos::null),
    xdot_(xdot),
    xdot_nc_(Teuchos::null),
    xdotdot_(xdotdot),
    xdotdot_nc_(Teuchos::null),
    stepperState_(stepperState),
    stepperState_nc_(Teuchos::null),
    physicsState_(physicsState),
    physicsState_nc_(Teuchos::null)
{
  if (stepperState_ == Teuchos::null) {
    stepperState_ = Teuchos::rcp(new StepperState<Scalar>("Default"));
  }
  if (physicsState_ == Teuchos::null) {
    physicsState_ = Teuchos::rcp(new PhysicsState<Scalar>());
  }
}

template <class Scalar>
SolutionState<Scalar>::SolutionState(const SolutionState<Scalar>& ss_)
  : metaData_(ss_.metaData_),
    metaData_nc_(ss_.metaData_nc_),
    x_(ss_.x_),
    x_nc_(ss_.x_nc_),
    xdot_(ss_.xdot_),
    xdot_nc_(ss_.xdot_nc_),
    xdotdot_(ss_.xdotdot_),
    xdotdot_nc_(ss_.xdotdot_nc_),
    stepperState_(ss_.stepperState_),
    stepperState_nc_(ss_.stepperState_nc_),
    physicsState_(ss_.physicsState_),
    physicsState_nc_(ss_.physicsState_nc_)
{
}

template <class Scalar>
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
  if (!Teuchos::is_null(stepperState_)) sS_out = stepperState_->clone();

  RCP<PhysicsState<Scalar> > pS_out;
  if (!Teuchos::is_null(physicsState_)) pS_out = physicsState_->clone();

  RCP<SolutionState<Scalar> > ss_out = Teuchos::rcp(new SolutionState<Scalar>(
      metaData_out, x_out, xdot_out, xdotdot_out, sS_out, pS_out));

  return ss_out;
}

template <class Scalar>
void SolutionState<Scalar>::copy(
    const Teuchos::RCP<const SolutionState<Scalar> >& ss)
{
  metaData_nc_->copy(ss->metaData_);
  this->copySolutionData(ss);
}

template <class Scalar>
void SolutionState<Scalar>::copySolutionData(
    const Teuchos::RCP<const SolutionState<Scalar> >& ss)
{
  if (ss->x_ == Teuchos::null)
    x_nc_ = Teuchos::null;
  else {
    if (x_nc_ == Teuchos::null) {
      x_nc_ = ss->x_->clone_v();
    }
    else
      Thyra::V_V(x_nc_.ptr(), *(ss->x_));
  }
  x_ = x_nc_;

  if (ss->xdot_ == Teuchos::null)
    xdot_nc_ = Teuchos::null;
  else {
    if (xdot_nc_ == Teuchos::null)
      xdot_nc_ = ss->xdot_->clone_v();
    else
      Thyra::V_V(xdot_nc_.ptr(), *(ss->xdot_));
  }
  xdot_ = xdot_nc_;

  if (ss->xdotdot_ == Teuchos::null)
    xdotdot_nc_ = Teuchos::null;
  else {
    if (xdotdot_nc_ == Teuchos::null)
      xdotdot_nc_ = ss->xdotdot_->clone_v();
    else
      Thyra::V_V(xdotdot_nc_.ptr(), *(ss->xdotdot_));
  }
  xdotdot_ = xdotdot_nc_;

  if (ss->stepperState_ == Teuchos::null)
    stepperState_nc_ = Teuchos::null;
  else {
    if (stepperState_nc_ == Teuchos::null)
      stepperState_nc_ = ss->stepperState_->clone();
    else
      stepperState_nc_->copy(ss->stepperState_);
  }
  stepperState_ = stepperState_nc_;

  if (ss->physicsState_ == Teuchos::null)
    physicsState_nc_ = Teuchos::null;
  else {
    if (physicsState_nc_ == Teuchos::null)
      physicsState_nc_ = ss->physicsState_->clone();
    else
      physicsState_nc_->copy(ss->physicsState_);
  }
  physicsState_ = physicsState_nc_;
}

template <class Scalar>
bool SolutionState<Scalar>::operator<(const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() < ss.metaData_->getTime());
}

template <class Scalar>
bool SolutionState<Scalar>::operator<=(const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() <= ss.metaData_->getTime());
}

template <class Scalar>
bool SolutionState<Scalar>::operator<(const Scalar& t) const
{
  return (this->metaData_->getTime() < t);
}

template <class Scalar>
bool SolutionState<Scalar>::operator<=(const Scalar& t) const
{
  return (this->metaData_->getTime() <= t);
}

template <class Scalar>
bool SolutionState<Scalar>::operator>(const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() > ss.metaData_->getTime());
}

template <class Scalar>
bool SolutionState<Scalar>::operator>=(const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() >= ss.metaData_->getTime());
}

template <class Scalar>
bool SolutionState<Scalar>::operator>(const Scalar& t) const
{
  return (this->metaData_->getTime() > t);
}

template <class Scalar>
bool SolutionState<Scalar>::operator>=(const Scalar& t) const
{
  return (this->metaData_->getTime() >= t);
}

template <class Scalar>
bool SolutionState<Scalar>::operator==(const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->getTime() == ss.metaData_->getTime());
}

template <class Scalar>
bool SolutionState<Scalar>::operator==(const Scalar& t) const
{
  return (this->metaData_->getTime() == t);
}

template <class Scalar>
std::string SolutionState<Scalar>::description() const
{
  std::ostringstream out;
  out << "SolutionState"
      << " (index =" << std::setw(6) << this->getIndex()
      << "; time =" << std::setw(10) << std::setprecision(3) << this->getTime()
      << "; dt =" << std::setw(10) << std::setprecision(3)
      << this->getTimeStep() << ")";
  return out.str();
}

template <class Scalar>
void SolutionState<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_EXTREME)) {
    metaData_->describe(*l_out, verbLevel);
    *l_out << "  x       = " << std::endl;
    x_->describe(*l_out, verbLevel);

    if (xdot_ != Teuchos::null) {
      *l_out << "  xdot_   = " << std::endl;
      xdot_->describe(*l_out, verbLevel);
    }
    if (xdotdot_ != Teuchos::null) {
      *l_out << "  xdotdot = " << std::endl;
      xdotdot_->describe(*l_out, verbLevel);
    }

    if (stepperState_ != Teuchos::null)
      stepperState_->describe(*l_out, verbLevel);
    if (physicsState_ != Teuchos::null)
      physicsState_->describe(*l_out, verbLevel);

    *l_out << std::string(this->description().length() + 8, '-') << std::endl;
  }
}

template <class Scalar>
void SolutionState<Scalar>::computeNorms(
    const Teuchos::RCP<const SolutionState<Scalar> >& ssIn)
{
  if (!getComputeNorms()) return;

  auto x = this->getX();
  this->setXNormL2(Thyra::norm(*x));

  if (ssIn != Teuchos::null) {
    auto xIn = ssIn->getX();

    // dx = x - xIn
    Teuchos::RCP<Thyra::VectorBase<Scalar> > dx =
        Thyra::createMember(x->space());
    Thyra::V_VmV(dx.ptr(), *x, *xIn);
    Scalar dxNorm  = Thyra::norm(*dx);
    Scalar xInNorm = Thyra::norm(*xIn);
    this->setDxNormL2Abs(dxNorm);
    // Compute change, e.g., ||x^n-x^(n-1)||/||x^(n-1)||
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar min = std::numeric_limits<Scalar>::min();
    if (xInNorm < min / eps) {  // numerically zero
      this->setDxNormL2Rel(std::numeric_limits<Scalar>::infinity());
    }
    else {
      // this->setDxNormL2Rel(dxNorm/(xInNorm + eps));
      this->setDxNormL2Rel(dxNorm / (xInNorm * (1.0 + 1.0e4 * eps)));
    }
  }
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > createSolutionStateX(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot)
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData_nc =
      Teuchos::rcp(new SolutionStateMetaData<Scalar>());

  Teuchos::RCP<StepperState<Scalar> > stepperState_nc =
      Teuchos::rcp(new StepperState<Scalar>("Default"));

  Teuchos::RCP<PhysicsState<Scalar> > physicsState_nc =
      Teuchos::rcp(new PhysicsState<Scalar>());

  Teuchos::RCP<SolutionState<Scalar> > ss = rcp(new SolutionState<Scalar>(
      metaData_nc, x, xdot, xdotdot, stepperState_nc, physicsState_nc));

  return ss;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > createSolutionStateX(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot)
{
  Teuchos::RCP<const SolutionStateMetaData<Scalar> > metaData =
      Teuchos::rcp(new SolutionStateMetaData<Scalar>());

  Teuchos::RCP<const StepperState<Scalar> > stepperState =
      Teuchos::rcp(new StepperState<Scalar>("Default"));

  Teuchos::RCP<const PhysicsState<Scalar> > physicsState =
      Teuchos::rcp(new PhysicsState<Scalar>());

  Teuchos::RCP<SolutionState<Scalar> > ss = rcp(new SolutionState<Scalar>(
      metaData, x, xdot, xdotdot, stepperState, physicsState));

  return ss;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > createSolutionStateME(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<StepperState<Scalar> >& stepperState,
    const Teuchos::RCP<PhysicsState<Scalar> >& physicsState)
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::rcp_const_cast;

  auto metaData_nc = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
  metaData_nc->setSolutionStatus(Status::PASSED);

  MEB::InArgs<Scalar> inArgs = model->getNominalValues();

  TEUCHOS_TEST_FOR_EXCEPTION(
      inArgs.supports(MEB::IN_ARG_x) == false, std::logic_error,
      model->description() << "does not support an x solution vector!");

  // The solution vector, x, is required (usually).
  auto x_nc = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());

  // The solution derivative, xdot, can be optional provided, based on
  // application needs.  Here we will base it on "supports" IN_ARG_x_dot.
  // Depending on the stepper used, a temporary xdot vector may be created
  // within the Stepper, but not moved to the SolutionState.
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot_nc;
  if (inArgs.supports(MEB::IN_ARG_x_dot)) {
    xdot_nc = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
  }
  else {
    xdot_nc = Teuchos::null;
  }

  // Similar as xdot.
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot_nc;
  if (inArgs.supports(MEB::IN_ARG_x_dot_dot)) {
    xdotdot_nc =
        rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot_dot());
  }
  else {
    xdotdot_nc = Teuchos::null;
  }

  Teuchos::RCP<StepperState<Scalar> > stepperState_nc;
  if (stepperState == Teuchos::null) {
    stepperState_nc = Teuchos::rcp(new StepperState<Scalar>());  // Use default
  }
  else {
    stepperState_nc = stepperState;
  }

  Teuchos::RCP<PhysicsState<Scalar> > physicsState_nc;
  if (physicsState == Teuchos::null) {
    physicsState_nc = Teuchos::rcp(new PhysicsState<Scalar>());  // Use default
  }
  else {
    physicsState_nc = physicsState;
  }

  Teuchos::RCP<SolutionState<Scalar> > ss =
      rcp(new SolutionState<Scalar>(metaData_nc, x_nc, xdot_nc, xdotdot_nc,
                                    stepperState_nc, physicsState_nc));

  return ss;
}

}  // namespace Tempus
#endif  // Tempus_SolutionState_impl_hpp
