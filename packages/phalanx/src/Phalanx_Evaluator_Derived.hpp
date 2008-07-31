
#ifndef PHX_FIELD_EVALUATOR_DERIVED_H
#define PHX_FIELD_EVALUATOR_DERIVED_H

#include <vector>

#include "Phalanx_Evaluator_Base.hpp"
#include "Phalanx_Evaluator_Utilities.hpp"

namespace PHX {

  template<typename EvalT, typename Traits>
  class EvaluatorDerived : 
    public PHX::EvaluatorBase<Traits> {
    
  public:
    
    EvaluatorDerived() {}

    virtual ~EvaluatorDerived() {}
    
  protected:
    
    PHX::EvaluatorUtilities<EvalT,Traits> utils;

  };

}

#endif
