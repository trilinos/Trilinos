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
#include "Panzer_Traits.hpp"

#include "Panzer_GatherTangent_Tpetra_decl.hpp"
#include "Panzer_GatherTangent_Tpetra_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherTangent_Tpetra,int,panzer::GlobalOrdinal)

#endif
