#ifndef TEMPUS_SOLUTIONSTATE_IMPL_HPP
#define TEMPUS_SOLUTIONSTATE_IMPL_HPP


namespace Tempus {

// SolutionState definitions:
template<class Scalar>
SolutionState<Scalar>::SolutionState()
{
  metaData_ = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState)
  : metaData_    (metaData),
    x_           (x),
    xdot_        (xdot),
    xdotdot_     (xdotdot),
    stepperState_(stepperState)
{}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Scalar time,
  const Scalar dt,
  const int    iStep,
  const Scalar errorAbs,
  const Scalar errorRel,
  const int    order,
  const int    nFailures,
  const int    nConsecutiveFailures,
  const Status solutionStatus,
  const bool   output,
  const bool   isInterpolated,
  const bool   isRestartable,
  const Scalar accuracy,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdot,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xdotdot,
  const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState)
  : x_           (x),
    xdot_        (xdot),
    xdotdot_     (xdotdot),
    stepperState_(stepperState)
{
  metaData_ = Teuchos::rcp(new SolutionStateMetaData<Scalar> (time,
      dt, iStep, errorAbs, errorRel, order, nFailures,
      nConsecutiveFailures, solutionStatus, output,
      isRestartable, isInterpolated, accuracy));
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(const SolutionState<Scalar>& ss_)
  :metaData_    (ss_.metaData_),
   x_           (ss_.x_),
   xdot_        (ss_.xdot_),
   xdotdot_     (ss_.xdotdot_),
   stepperState_(ss_.stepperState_)
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

  RCP<StepperState<Scalar> > stepperState_out;
  if (!Teuchos::is_null(stepperState_)) stepperState_out=stepperState_->clone();

  RCP<SolutionState<Scalar> > ss_out = Teuchos::rcp(new SolutionState<Scalar> (
    metaData_out, x_out, xdot_out, xdotdot_out, stepperState_out));

  return ss_out;
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->time < ss->metaData_->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->time <= ss.metaData_->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const Scalar& t) const
{
  return (this->metaData_->time < t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const Scalar& t) const
{
  return (this->metaData_->time <= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->time > ss.metaData_->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->time >= ss.metaData_->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const Scalar& t) const
{
  return (this->metaData_->time > t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const Scalar& t) const
{
  return (this->metaData_->time >= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const SolutionState<Scalar>& ss) const
{
  return (this->metaData_->time == ss.metaData_->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const Scalar& t) const
{
  return (this->metaData_->time == t);
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
  }
}

} // namespace Tempus
#endif // TEMPUS_SOLUTIONSTATE_IMPL_HPP
