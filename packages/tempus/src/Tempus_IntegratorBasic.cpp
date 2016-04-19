#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorBasic_impl.hpp"

namespace Tempus {

template class IntegratorBasic<double>;

template Teuchos::RCP<IntegratorBasic<double> >
IntegratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>     parameterList,
  const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
  const Teuchos::RCP<Thyra::VectorBase<double> >&x,
  const Teuchos::RCP<Thyra::VectorBase<double> >&xdot=Teuchos::null,
  const Teuchos::RCP<Thyra::VectorBase<double> >&xdotdot=Teuchos::null);

} // namespace Tempus
