
#ifndef PHX_EXAMPLE_VP_DENSITY_HPP
#define PHX_EXAMPLE_VP_DENSITY_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename EvalT, typename Traits>
class Density : 
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  Density(const Teuchos::ParameterList& p);
  
  ~Density();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData ud);
  
  void preEvaluate() {}
  
  void postEvaluate() {}
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  double constant;

  PHX::Field<ScalarT> density;
  PHX::Field<ScalarT> temp;

  std::size_t data_layout_size;

};

#include "Evaluator_Density_Def.hpp"

#endif
