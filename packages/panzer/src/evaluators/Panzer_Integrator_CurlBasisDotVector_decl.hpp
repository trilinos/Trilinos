#ifndef PANZER_EVALUATOR_CURLBASISDOTVECTOR_DECL_HPP
#define PANZER_EVALUATOR_GRADBASISDOTVECTOR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
/** In 3D this computes
  * 
  *  \f$\int \nabla\times \phi \cdot v \f$
  *
  * however the name can be misleading. The curl of a vector
  * in 2D is simply a scalar, here the evaluators handles
  * both cases.
  */
PHX_EVALUATOR_CLASS(Integrator_CurlBasisDotVector)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;
    
  PHX::MDField<ScalarT> flux; // note that this is used for both vector and scalar fields
                              // but this follows the distinction between 2D and 3D curls
    
  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;
  std::size_t num_qp;
  std::size_t num_dim;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  bool useScalarField;

  Intrepid::FieldContainer<ScalarT> tmp;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
PHX_EVALUATOR_CLASS_END

}

#endif
