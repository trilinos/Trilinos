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

#include "Panzer_ScatterResidual_BlockedEpetra_decl.hpp"
#include "Panzer_ScatterResidual_BlockedEpetra_impl.hpp"

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
#include "Panzer_ScatterResidual_BlockedEpetra_Hessian_impl.hpp"
#endif

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterResidual_BlockedEpetra,int,int)