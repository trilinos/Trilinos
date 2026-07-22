// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_ALIAS_FIELD_HPP
#define PHX_EVALUATOR_ALIAS_FIELD_HPP

#include "Phalanx_config.hpp"

#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif

namespace PHX {

  //! Evaluator to help set dependencies for aliased fields
  template<typename EvalT, typename Traits>
  class AliasField :
#ifdef PHX_ENABLE_KOKKOS_AMT
    public PHX::TaskBase<Traits,PHX::AliasField<EvalT,Traits> >
#else
    public PHX::EvaluatorWithBaseImpl<Traits>
#endif
  {
    
  public:
    AliasField(const PHX::FieldTag& aliasedField,
               const PHX::FieldTag& targetField)
    {
      this->addEvaluatedField(aliasedField);
      this->addDependentField(targetField);
      this->setName("Alias Field: \"" + aliasedField.name()
                    + "\" to \"" + targetField.name() + "\"");
    }
    
    void postRegistrationSetup(typename Traits::SetupData ,
                               PHX::FieldManager<Traits>& ) {}
    
    void evaluateFields(typename Traits::EvalData ) {}

    KOKKOS_INLINE_FUNCTION
    void operator () (const int ) const {}
  };
  
}

#endif
