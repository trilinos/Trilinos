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

#include "Panzer_ExplicitModelEvaluator.hpp"
#include "Panzer_ExplicitModelEvaluator_impl.hpp"

namespace panzer {

template class ExplicitModelEvaluator<panzer::Traits::RealType>; 

}
