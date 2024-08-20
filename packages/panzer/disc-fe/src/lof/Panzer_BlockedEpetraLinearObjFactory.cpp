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

#include "Panzer_Traits.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory_impl.hpp"

template class panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>;