// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_IntegratorForwardSensitivity_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorForwardSensitivity)

  // non-member ctor
  template Teuchos::RCP<IntegratorForwardSensitivity<double> >
  integratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList>        parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

  // non-member ctor
  template Teuchos::RCP<IntegratorForwardSensitivity<double> >
  integratorForwardSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    std::string stepperType);

  // non-member ctor
  template Teuchos::RCP<IntegratorForwardSensitivity<double> >
  integratorForwardSensitivity();

} // namespace Tempus

#endif
