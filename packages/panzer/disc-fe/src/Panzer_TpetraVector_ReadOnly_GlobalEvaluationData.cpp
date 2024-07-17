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

#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData_impl.hpp"

#include "Panzer_NodeType.hpp"

// template class panzer::TpetraVector_ReadOnly_GlobalEvaluationData<double,int,int>;
template class panzer::TpetraVector_ReadOnly_GlobalEvaluationData<double,int,panzer::GlobalOrdinal>;
