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
#include "Panzer_ExplicitTemplateInstantiation.hpp"

// Files for this specific example.
#include "mySourceTerm.hpp"
#include "mySourceTermImpl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(MySourceTerm)

// end of mySourceTerm.cpp
