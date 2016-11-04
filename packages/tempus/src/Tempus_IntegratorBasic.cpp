#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorBasic_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorBasic)

  // non-member ctor
  template Teuchos::RCP<IntegratorBasic<double> >
  integratorBasic(Teuchos::RCP<Teuchos::ParameterList>        parameterList,
                  const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

} // namespace Tempus

#endif
