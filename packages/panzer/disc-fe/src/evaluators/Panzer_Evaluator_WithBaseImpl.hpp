// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_WITHBASEIMPL_HPP
#define PANZER_EVALUATOR_WITHBASEIMPL_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Panzer_Evaluator_DomainInterface.hpp"
#include "Panzer_Workset.hpp"

namespace panzer {

//! Wrapper to PHX::EvaluatorWithBaseImpl that implements Panzer-specific helpers.
template<typename TRAITS>
class EvaluatorWithBaseImpl :
    public PHX::EvaluatorWithBaseImpl<TRAITS>,
    public panzer::DomainEvaluator {

public:
  //! An evaluator builder sets the details index.
  void setDetailsIndex(const int di) { wda.setDetailsIndex(di); }
  
protected:
  WorksetDetailsAccessor wda;
};

}

#endif
