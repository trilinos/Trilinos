
#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Utilities.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename ScalarT, typename Traits>
class NonlinearSource : public PHX::EvaluatorUtilities<Traits>,
			public PHX::EvaluatorDerived<ScalarT, Traits> {
  
public:
  
  NonlinearSource(const Teuchos::ParameterList& p);
  
  ~NonlinearSource();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
  void preEvaluate(typename Traits::PreEvalData d);
  
  void postEvaluate(typename Traits::PostEvalData d);
  
private:
  
  PHX::Field<ScalarT> source;
  PHX::Field<ScalarT> density;
  PHX::Field<ScalarT> temp;

  std::size_t data_layout_size;
};

#include "Evaluator_NonlinearSource_Def.hpp"

#endif
