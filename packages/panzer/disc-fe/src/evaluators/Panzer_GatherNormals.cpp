// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_ExplicitTemplateInstantiation.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_GatherNormals_decl.hpp"
#include "Panzer_GatherNormals_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::GatherNormals)
