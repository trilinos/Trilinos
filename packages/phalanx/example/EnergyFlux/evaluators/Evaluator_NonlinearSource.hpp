
#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename EvalT, typename Traits>
class NonlinearSource : public PHX::EvaluatorWithBaseImpl<Traits>,
			public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  NonlinearSource(const Teuchos::ParameterList& p);
  
  ~NonlinearSource();

  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
  void preEvaluate(typename Traits::PreEvalData d);
  
  void postEvaluate(typename Traits::PostEvalData d);
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  PHX::Field<ScalarT> source;
  PHX::Field<ScalarT> density;
  PHX::Field<ScalarT> temp;

  std::size_t data_layout_size;
};

#include "Evaluator_NonlinearSource_Def.hpp"

#endif
