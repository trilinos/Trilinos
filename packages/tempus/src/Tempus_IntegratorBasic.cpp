#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorBasic_impl.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION

namespace Tempus {

template class IntegratorBasic<double>;

template Teuchos::RCP<IntegratorBasic<double> >
IntegratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                     parameterList,
  const Teuchos::RCP<Thyra::ModelEvaluator<double> >&      model,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<double> >& solver);

} // namespace Tempus

#endif
