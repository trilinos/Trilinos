// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "PanzerDiscFE_config.hpp"
#include "Panzer_ExplicitTemplateInstantiation.hpp"
#include "Panzer_ZeroContributedField.hpp"
#include "Panzer_ZeroContributedField_impl.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Instantiate the class.
//
///////////////////////////////////////////////////////////////////////////////
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::ZeroContributedField)
