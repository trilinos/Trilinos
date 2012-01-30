#ifndef PANZER_EVALUATOR_GRADBASISDOTVECTOR_DECL_HPP
#define PANZER_EVALUATOR_GRADBASISDOTVECTOR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS(Integrator_GradBasisDotVector)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
    
  PHX::MDField<ScalarT,Cell,IP,Dim> flux;
    
  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;

  std::size_t num_qp;

  std::size_t num_dim;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  Intrepid::FieldContainer<ScalarT> tmp;

PHX_EVALUATOR_CLASS_END

}

#endif
