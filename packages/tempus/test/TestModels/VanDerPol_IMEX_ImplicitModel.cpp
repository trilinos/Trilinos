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
#include "VanDerPol_IMEX_ImplicitModel.hpp"
#include "VanDerPol_IMEX_ImplicitModel_impl.hpp"

namespace Tempus_Test {
TEMPUS_INSTANTIATE_TEMPLATE_CLASS(VanDerPol_IMEX_ImplicitModel)
}

#endif
