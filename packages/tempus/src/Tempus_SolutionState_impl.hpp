#ifndef TEMPUS_SOLUTIONSTATE_IMPL_HPP
#define TEMPUS_SOLUTIONSTATE_IMPL_HPP

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

namespace Tempus {

// SolutionState definitions:
template<class Scalar>
SolutionState<Scalar>::SolutionState()
{
  metaData = Teuchos::rcp(new SolutionStateMetaData<Scalar>());
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const RCP<SolutionStateMetaData<Scalar> > metaData_,
  const RCP<Thyra::VectorBase<Scalar> >& x_,
  const RCP<Thyra::VectorBase<Scalar> >& xdot_,
  const RCP<Thyra::VectorBase<Scalar> >& xdotdot_,
  const RCP<Tempus::StepperState<Scalar> >& stepperState_)
  : metaData     (metaData_),
    x            (x_),
    xdot         (xdot_),
    xdotdot      (xdotdot_),
    stepperState (stepperState_)
{}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Scalar time_,
  const Scalar dt_,
  const int    iStep_,
  const Scalar errorAbs_,
  const Scalar errorRel_,
  const int    order_,
  const int    nFailures_,
  const int    nConsecutiveFailures_,
  const SolutionStatus status_,
  const bool   output_,
  const bool   isAccepted_,
  const bool   isInterpolated_,
  const bool   isRestartable_,
  const Scalar accuracy_,
  const RCP<Thyra::VectorBase<Scalar> >& x_,
  const RCP<Thyra::VectorBase<Scalar> >& xdot_,
  const RCP<Thyra::VectorBase<Scalar> >& xdotdot_,
  const RCP<Tempus::StepperState<Scalar> >& stepperState_)
  : x            (x_),
    xdot         (xdot_),
    xdotdot      (xdotdot_),
    stepperState (stepperState_)
{
  metaData = Teuchos::rcp(new SolutionStateMetaData<Scalar> (time_,
      dt_, iStep_, errorAbs_, errorRel_, order_, nFailures_,
      nConsecutiveFailures_, status_, output_, isAccepted_,
      isRestartable_, isInterpolated_, accuracy_));
}

template<class Scalar>
SolutionState<Scalar>::SolutionState(const SolutionState<Scalar>& ss_)
  :metaData     (ss_.metaData),
   x            (ss_.x),
   xdot         (ss_.xdot),
   xdotdot      (ss_.xdotdot),
   stepperState (ss_.stepperState)
{}

template<class Scalar>
RCP<SolutionState<Scalar> > SolutionState<Scalar>::clone() const
{
  RCP<Thyra::VectorBase<Scalar> > x_out;
  if (!Teuchos::is_null(x)) x_out = x->clone_v();

  RCP<Thyra::VectorBase<Scalar> > xdot_out;
  if (!Teuchos::is_null(xdot)) xdot_out = xdot->clone_v();

  RCP<Thyra::VectorBase<Scalar> > xdotdot_out;
  if (!Teuchos::is_null(xdotdot)) xdotdot_out = xdotdot->clone_v();

  RCP<StepperState<Scalar> > stepperState_out;
  if (!Teuchos::is_null(stepperState)) stepperState_out=stepperState->clone();

  RCP<SolutionState<Scalar> > ss_out = Teuchos::rcp(new SolutionState<Scalar> (
    metaData, x_out, xdot_out, xdotdot_out, stepperState_out));

  return ss_out;
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const SolutionState<Scalar>& ss) const
{
  return (this->metaData->time < ss->metaData->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData->time <= ss.metaData->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const Scalar& t) const
{
  return (this->metaData->time < t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const Scalar& t) const
{
  return (this->metaData->time <= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const SolutionState<Scalar>& ss) const
{
  return (this->metaData->time > ss.metaData->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const SolutionState<Scalar>& ss) const
{
  return (this->metaData->time >= ss.metaData->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const Scalar& t) const
{
  return (this->metaData->time > t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const Scalar& t) const
{
  return (this->metaData->time >= t);
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const SolutionState<Scalar>& ss) const
{
  return (this->metaData->time == ss.metaData->time);
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const Scalar& t) const
{
  return (this->metaData->time == t);
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
        metaData->describe(out,verbLevel);
    out << "x = " << std::endl;
    x->describe(out,verbLevel);
    if (xdot != Teuchos::null) {
      out << "xdot = " << std::endl;
      xdot->describe(out,verbLevel);
    }
    if (xdotdot != Teuchos::null) {
      out << "xdotdot = " << std::endl;
      xdotdot->describe(out,verbLevel);
    }
    if (stepperState != Teuchos::null) {
      out << "stepperState = " << std::endl;
      xdotdot->describe(out,verbLevel);
    }
  }
}

} // namespace Tempus
#endif // TEMPUS_SOLUTIONSTATE_IMPL_HPP
