// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerDiscFE_config.hpp"

#include "Panzer_ExplicitTemplateInstantiation.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_ScatterDirichletResidual_BlockedTpetra.hpp"
#include "Panzer_ScatterDirichletResidual_BlockedTpetra_impl.hpp"

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterDirichletResidual_BlockedTpetra_Hessian_impl.hpp"
#endif

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterDirichletResidual_BlockedTpetra,int,panzer::GlobalOrdinal)
