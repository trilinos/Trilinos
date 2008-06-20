
#ifndef PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP
#define PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldEvaluator_Utilities.hpp"
#include "Phalanx_FieldEvaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename ScalarT, typename Traits>
class FEInterpolation : public PHX::FieldEvaluatorUtilities<Traits>,
			public PHX::FieldEvaluatorDerived<ScalarT, Traits>  {
  
public:
  
  FEInterpolation(const Teuchos::ParameterList& p);
  
  ~FEInterpolation();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
  void preEvaluate() {}
  
  void postEvaluate() {}
  
private:

  std::vector< PHX::Field< ScalarT> > s_n;
  std::vector< PHX::Field< ScalarT> > s_qp;
  std::vector< PHX::Field< ScalarT> > g_n;
  std::vector< PHX::Field< MyVector<ScalarT> > > g_qp;
  
};

#include "Evaluator_FEInterpolation_Def.hpp"

#endif
