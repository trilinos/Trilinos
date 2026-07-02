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
  template Teuchos::RCP<StepperExponential_EPI2<double> >
  createStepperExponential_EPI2(
      const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
      Teuchos::RCP<Teuchos::ParameterList> pl);

  template Teuchos::RCP<StepperExponential_EPI3<double> >
  createStepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl);

}

#endif
