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

#include "Panzer_ResponseScatterEvaluator_Probe.hpp"
#include "Panzer_ResponseScatterEvaluator_Probe_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ResponseScatterEvaluator_Probe,int,panzer::GlobalOrdinal)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ResponseScatterEvaluator_ProbeBase,int,panzer::GlobalOrdinal)
