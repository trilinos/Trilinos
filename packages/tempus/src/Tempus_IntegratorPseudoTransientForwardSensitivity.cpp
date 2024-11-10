//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorPseudoTransientForwardSensitivity)

/// Nonmember constructor
template Teuchos::RCP<IntegratorPseudoTransientForwardSensitivity<double> >
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& sens_solve_model);

/// Nonmember constructor
template Teuchos::RCP<IntegratorPseudoTransientForwardSensitivity<double> >
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& sens_residual_model);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientForwardSensitivity<double> >
createIntegratorPseudoTransientForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientForwardSensitivity<double> >
createIntegratorPseudoTransientForwardSensitivity();

}  // namespace Tempus

#endif
