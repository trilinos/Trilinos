// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Traits.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

panzer::Traits::PED::PED():
  gedc(new GlobalEvaluationDataContainer())
{
  // This only exists to initialize GlobalEvaluationDataContainer
}
