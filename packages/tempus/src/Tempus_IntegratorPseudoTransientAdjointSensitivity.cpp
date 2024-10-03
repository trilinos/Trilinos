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
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorPseudoTransientAdjointSensitivity)

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    std::string stepperType);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_model);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_model,
    std::string stepperType);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_solve_model);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_solve_model,
    std::string stepperType);

// Nonmember ctor
template Teuchos::RCP<IntegratorPseudoTransientAdjointSensitivity<double> >
integratorPseudoTransientAdjointSensitivity();

}  // namespace Tempus

#endif
