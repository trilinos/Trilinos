
#ifndef PHX_EXAMPLE_VP_FOURIER_HPP
#define PHX_EXAMPLE_VP_FOURIER_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "AlgebraicTypes.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DataLayout_Generic.hpp"
#include "Phalanx_Field.hpp"

template<typename EvalT, typename Traits>
class Fourier : public PHX::EvaluatorWithBaseImpl<Traits>,
		public PHX::EvaluatorDerived<EvalT, Traits>  {
  
public:
  
  Fourier(const Teuchos::ParameterList& p);
  
  ~Fourier();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
  void preEvaluate() {}
  
  void postEvaluate() {}
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  PHX::Field< MyVector<ScalarT> > flux;
  PHX::Field< ScalarT > density;
  PHX::Field< ScalarT > dc;
  PHX::Field< MyVector<ScalarT> > grad_temp;
  
  std::size_t data_layout_size;

};

#include "Evaluator_EnergyFlux_Fourier_Def.hpp"

#endif
