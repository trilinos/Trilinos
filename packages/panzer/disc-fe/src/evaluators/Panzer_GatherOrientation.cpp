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

#include "Panzer_GatherOrientation_decl.hpp"
#include "Panzer_GatherOrientation_impl.hpp"

// This is limited to only single value scalar types, the blocked
// version must use the vector input
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherOrientation,int,panzer::GlobalOrdinal)
