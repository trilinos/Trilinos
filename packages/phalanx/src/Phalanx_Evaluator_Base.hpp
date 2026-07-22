// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_EVALUATOR_BASE_HPP
#define PHX_FIELD_EVALUATOR_BASE_HPP

namespace PHX {
  

  /*! \brief Template Manager "Base" class object for all field evaluators.
  */
  template<typename Traits>
  class EvaluatorBase {
    
  public:
    
    EvaluatorBase() {}
    
    virtual ~EvaluatorBase() {}
    
  };
  
}

#endif 
