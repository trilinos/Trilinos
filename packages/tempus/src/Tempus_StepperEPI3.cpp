// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_StepperEPI3.hpp"
#include "Tempus_StepperEPI3_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperEPI3)

  // Nonmember constructor
  template Teuchos::RCP<StepperEPI3<double> >
  createStepperEPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}

#endif
