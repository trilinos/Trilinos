
#ifndef PHX_EXAMPLE_VP_CONSTANT_HPP
#define PHX_EXAMPLE_VP_CONSTANT_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Utilities.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename ScalarT, typename Traits>
class Constant : 
  public PHX::EvaluatorUtilities<Traits>,
  public PHX::EvaluatorDerived<ScalarT, Traits> {
  
public:
  
  Constant(Teuchos::ParameterList& p);
  
  ~Constant();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData ud);
  
  void preEvaluate() {}
  
  void postEvaluate() {}
  
private:
  
  ScalarT value;

  PHX::Field<ScalarT> constant;

};

#include "Evaluator_Constant_Def.hpp"

#endif
