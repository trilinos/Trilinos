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
#include "Tempus_IntegratorAdjointSensitivity.hpp"
#include "Tempus_IntegratorAdjointSensitivity_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorAdjointSensitivity)

// Nonmember ctor
template Teuchos::RCP<IntegratorAdjointSensitivity<double> >
createIntegratorAdjointSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& adjoint_model);

// Nonmember ctor
template Teuchos::RCP<IntegratorAdjointSensitivity<double> >
createIntegratorAdjointSensitivity();

}  // namespace Tempus

#endif
