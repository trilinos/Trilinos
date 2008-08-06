// @HEADER
// @HEADER

#ifndef PHX_EXAMPLE_VP_CONSTANT_HPP
#define PHX_EXAMPLE_VP_CONSTANT_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

namespace Teuchos {
  class ParameterList;
}

template<typename EvalT, typename Traits>
class Constant : 
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  Constant(Teuchos::ParameterList& p);
  
  ~Constant();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData ud);
  
  void preEvaluate() {}
  
  void postEvaluate() {}
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  ScalarT value;

  PHX::Field<ScalarT> constant;

};

#include "Evaluator_Constant_Def.hpp"

#endif
