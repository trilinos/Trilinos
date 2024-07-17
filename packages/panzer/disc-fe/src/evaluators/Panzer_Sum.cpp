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

#include "Panzer_Sum.hpp"
#include "Panzer_Sum_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::Sum)

PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(panzer::SumStatic,panzer::Cell)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::SumStatic,panzer::Cell,panzer::IP)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::SumStatic,panzer::Cell,panzer::BASIS)
