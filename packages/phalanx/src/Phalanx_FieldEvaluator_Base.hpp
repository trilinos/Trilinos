#ifndef PHX_FIELD_EVALUATOR_BASE_HPP
#define PHX_FIELD_EVALUATOR_BASE_HPP

namespace PHX {
  

  /*! \brief Template Manager "Base" class object for all field evaluators.
  */
  template<typename Traits>
  class FieldEvaluatorBase {
    
  public:
    
    FieldEvaluatorBase() {}
    
    virtual ~FieldEvaluatorBase() {}
    
  };
  
}

#endif 
