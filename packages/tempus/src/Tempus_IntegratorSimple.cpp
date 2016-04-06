#include "Tempus_IntegratorSimple.hpp"
#include "Tempus_IntegratorSimple_impl.hpp"

namespace tempus {

template class IntegratorSimple< double >;

template Teuchos::RCP< IntegratorSimple< double > >
IntegratorSimple(
  RCP<ParameterList>              parameterList,
  RCP<Thyra::VectorBase<double> > x,
  RCP<Thyra::VectorBase<double> > xdot=Teuchos::null,
  RCP<Thyra::VectorBase<double> > xdotdot=Teuchos::null );

} // namespace tempus
