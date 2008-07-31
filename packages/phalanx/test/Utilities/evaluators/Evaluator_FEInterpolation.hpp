
#ifndef PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP
#define PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Field.hpp"

/** \brief Finite Element Interpolation Evaluator

    This object evaluates a scalar field and it's gradient at the
    quadrature points for a specific variable.

*/
template<typename EvalT, typename Traits>
class FEInterpolation : public PHX::EvaluatorWithBaseImpl<Traits>,
			public PHX::EvaluatorDerived<EvalT, Traits>  {
  
public:
  
  FEInterpolation(const Teuchos::ParameterList& p);
  
  ~FEInterpolation();
  
  void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  //! Values at nodes
  PHX::Field< ScalarT > val_node;

  //! Values at quadrature points
  PHX::Field< ScalarT > val_qp;

  //! Gradient values at quadrature points
  PHX::Field< MyVector<ScalarT> > val_grad_qp;
  
};

#include "Evaluator_FEInterpolation_Def.hpp"

#endif
