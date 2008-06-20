
#ifndef PHX_FIELD_EVALUATOR_DERIVED_H
#define PHX_FIELD_EVALUATOR_DERIVED_H

#include <vector>

#include "Phalanx_FieldEvaluator_Base.hpp"

namespace PHX {

  template<typename ScalarT, typename Traits>
  class FieldEvaluatorDerived : 
    public PHX::FieldEvaluatorBase<Traits> {
    
  public:
    
    FieldEvaluatorDerived() {}

    virtual ~FieldEvaluatorDerived() {}
    
  };

}

#endif
