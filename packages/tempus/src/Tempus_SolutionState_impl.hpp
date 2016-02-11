#ifndef TEMPUS_SOLUTIONSTATE_IMPL_HPP
#define TEMPUS_SOLUTIONSTATE_IMPL_HPP

#include "Tempus_SolutionState.hpp"

namespace tempus {

// SolutionState definitions:
template<class Scalar>
SolutionState<Scalar>::SolutionState()
  :time(-1),
   dt(-1),
   dtMin(-1),
   dtMax(-1),
   iStep(-1),
   order(-1),
   error(-1),
   isInterpolated(false),
   accuracy(-1)
{}

template<class Scalar>
SolutionState<Scalar>::SolutionState(
  const Scalar time_,
  const Scalar dt_,
  const Scalar dtMin_,
  const Scalar dtMax_,
  const int    iStep_,
  const int    order_,
  const Scalar error_,
  const bool   isInterpolated_,
  const Scalar accuracy_,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& x_,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdot_,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& xdotdot_)
  :time          (time_),
   dt            (dt_),
   dtMin         (dtMin_),
   dtMax         (dtMax_),
   iStep         (iStep_),
   order         (order_),
   error         (error_),
   isInterpolated(isInterpolated_),
   accuracy      (accuracy_),
   x             (x_),
   xdot          (xdot_),
   xdotdot       (xdotdot_),
{}

template<class Scalar>
SolutionState<Scalar>::SolutionState( const SolutionState<Scalar>& ss_ )
  :time          (ss_.time),
   dt            (ss_.dt),
   dtMin         (ss_.dtMin),
   dtMax         (ss_.dtMax),
   iStep         (ss_.iStep),
   order         (ss_.order),
   error         (ss_.error),
   isInterpolated(ss_.isInterpolated_),
   accuracy      (ss_.accuracy),
   x             (ss_.x),
   xdot          (ss_.xdot),
   xdotdot       (ss_.xdotdot),
{}

template<class Scalar>
RCP<SolutionState<Scalar> > SolutionState<Scalar>::clone() const
{
  RCP<VectorBase<Scalar> > x_out;
  if (!Teuchos::is_null(x)) x_out = x->clone_v();

  RCP<VectorBase<Scalar> > xdot_out;
  if (!Teuchos::is_null(xdot)) xdot_out = xdot->clone_v();

  RCP<VectorBase<Scalar> > xdotdot_out;
  if (!Teuchos::is_null(xdotdot)) xdotdot_out = xdotdot->clone_v();

  RCP<SolutionState<Scalar> > ss_out = Teuchos::rcp(new SolutionState<Scalar> (
    time, dt, dtMin, dtMax, iStep, order, error, isInterpolated, accuracy,
    x_out, xdot_out, xdotdot_out));

  return ss_out;
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const SolutionState<Scalar>& ss) const
{
  return( this->time < ss.time );
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const SolutionState<Scalar>& ss) const
{
  return( this->time <= ss.time );
}

template<class Scalar>
bool SolutionState<Scalar>::operator< (const Scalar& t) const
{
  return( this->time < t );
}

template<class Scalar>
bool SolutionState<Scalar>::operator<= (const Scalar& t) const
{
  return( this->time <= t );
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const SolutionState<Scalar>& ss) const
{
  return( this->time > ss.time );
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const SolutionState<Scalar>& ss) const
{
  return( this->time >= ss.time );
}

template<class Scalar>
bool SolutionState<Scalar>::operator> (const Scalar& t) const
{
  return( this->time > t );
}

template<class Scalar>
bool SolutionState<Scalar>::operator>= (const Scalar& t) const
{
  return( this->time >= t );
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const SolutionState<Scalar>& ss) const
{
  return( this->time == ss.time );
}

template<class Scalar>
bool SolutionState<Scalar>::operator== (const Scalar& t) const
{
  return( this->time == t );
}

template<class Scalar>
std::string SolutionState<Scalar>::description() const
{
  std::string name = "Tempus::SolutionState";
  return(name);
}

template<class Scalar>
void SolutionState<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  if (verbLevel == Teuchos::VERB_EXTREME) {
    out << description() << "::describe:" << std::endl;
    out << "time           = " << time << std:endl;
    out << "dt             = " << dt << std:endl;
    out << "dtMin          = " << dtMin << std:endl;
    out << "dtMax          = " << dtMax << std:endl;
    out << "iStep          = " << iStep << std:endl;
    out << "order          = " << order << std:endl;
    out << "error          = " << error << std:endl;
    out << "isInterpolated = " << isInterpolated << std:endl;
    out << "accuracy       = " << accuracy << std:endl;
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
  }
}

} // namespace tempus
#endif TEMPUS_SOLUTIONSTATE_IMPL_HPP
