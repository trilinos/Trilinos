#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorBasic_impl.hpp"

namespace tempus {

template class IntegratorBasic< double >;

template Teuchos::RCP< IntegratorBasic< double > >
IntegratorBasic(
  RCP<ParameterList>              parameterList,
  RCP<Thyra::VectorBase<double> > x,
  RCP<Thyra::VectorBase<double> > xdot=Teuchos::null,
  RCP<Thyra::VectorBase<double> > xdotdot=Teuchos::null );

} // namespace tempus
