#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionState_impl.hpp"

namespace tempus {

template class SolutionState< double >;

template Teuchos::RCP< SolutionState< double > > SolutionState();

template Teuchos::RCP< SolutionState< double > >
SolutionState( const double time,
               const double dt,
               const double dtMin,
               const double dtMax,
               const int    iStep,
               const int    order,
               const double error,
               const bool   isInterpolated,
               const bool   isRestartable,
               const double accuracy,
               const Teuchos::RCP<const Thyra::VectorBase<double> >& x,
               const Teuchos::RCP<const Thyra::VectorBase<double> >& xdot,
               const Teuchos::RCP<const Thyra::VectorBase<double> >& xdotdot);

} // namespace tempus
