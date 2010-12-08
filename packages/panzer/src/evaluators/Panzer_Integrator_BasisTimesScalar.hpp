#ifndef PANZER_EVALUATOR_BASISTIMESSCALAR_HPP
#define PANZER_EVALUATOR_BASISTIMESSCALAR_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS(Integrator_BasisTimesScalar)
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
    
  PHX::MDField<ScalarT,Cell,IP> scalar;

  std::size_t num_nodes;

  std::size_t num_qp;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  Intrepid::FieldContainer<ScalarT> tmp;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_Integrator_BasisTimesScalarT.hpp"

#endif
