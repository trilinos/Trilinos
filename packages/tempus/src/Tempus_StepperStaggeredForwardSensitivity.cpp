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
#include "Tempus_StepperStaggeredForwardSensitivity_decl.hpp"
#include "Tempus_StepperStaggeredForwardSensitivity_impl.hpp"

namespace Tempus {
TEMPUS_INSTANTIATE_TEMPLATE_CLASS(StepperStaggeredForwardSensitivity)
}

#endif
