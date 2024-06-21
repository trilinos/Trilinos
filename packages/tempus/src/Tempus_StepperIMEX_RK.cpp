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
#include "Tempus_StepperIMEX_RK.hpp"
#include "Tempus_StepperIMEX_RK_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperIMEX_RK)

// Nonmember constructor
template Teuchos::RCP<StepperIMEX_RK<double> > createStepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> >& model,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif
