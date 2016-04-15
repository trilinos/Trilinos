#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionState_impl.hpp"

namespace tempus {

template class SolutionState<double>;

template RCP<SolutionState<double> > SolutionState();

template RCP<SolutionState<double> >
SolutionState(const RCP<SolutionStateMetaData<double> > ssmd,
              const RCP<Thyra::VectorBase<double> >& x,
              const RCP<Thyra::VectorBase<double> >& xdot,
              const RCP<Thyra::VectorBase<double> >& xdotdot,
              const RCP<tempus::StepperState<double> >& stepperState);

template RCP<SolutionState<double> >
SolutionState(const double time,
              const double dt,
              const int    iStep,
              const double errorAbs,
              const double errorRel,
              const int    order,
              const int    nFailures,
              const int    nConsecutiveFailures,
              const SolutionStatus status,
              const bool   output,
              const bool   isAccepted,
              const bool   isRestartable,
              const bool   isInterpolated,
              const double accuracy,
              const RCP<Thyra::VectorBase<double> >& x,
              const RCP<Thyra::VectorBase<double> >& xdot,
              const RCP<Thyra::VectorBase<double> >& xdotdot,
              const RCP<tempus::StepperState<double> >& stepperState);

} // namespace tempus
