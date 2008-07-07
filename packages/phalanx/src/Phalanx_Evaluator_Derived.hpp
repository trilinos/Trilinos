
#ifndef PHX_FIELD_EVALUATOR_DERIVED_H
#define PHX_FIELD_EVALUATOR_DERIVED_H

#include <vector>

#include "Phalanx_Evaluator_Base.hpp"

namespace PHX {

  template<typename ScalarT, typename Traits>
  class EvaluatorDerived : 
    public PHX::EvaluatorBase<Traits> {
    
  public:
    
    EvaluatorDerived() {}

    virtual ~EvaluatorDerived() {}
    
  };

}

#endif
