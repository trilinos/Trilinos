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

  // Nonmember ctor
  template Teuchos::RCP<IntegratorForwardSensitivity<double> >
  createIntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList>        parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

  // Nonmember ctor
  template Teuchos::RCP<IntegratorForwardSensitivity<double> >
  createIntegratorForwardSensitivity();

} // namespace Tempus

#endif
