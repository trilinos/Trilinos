// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_StepperEPI.hpp"
#include "Tempus_StepperEPI_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperEPI)

  // Nonmember constructor
  template Teuchos::RCP<StepperEPI<double> >
  createStepperEPI(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}

#endif
