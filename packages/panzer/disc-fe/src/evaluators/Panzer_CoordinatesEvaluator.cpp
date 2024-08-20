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

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_CoordinatesEvaluator.hpp"
#include "Panzer_CoordinatesEvaluator_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::CoordinatesEvaluator)

#endif
