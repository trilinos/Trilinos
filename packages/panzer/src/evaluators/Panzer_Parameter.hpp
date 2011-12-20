#ifndef PANZER_EVALUATOR_PARAMETER_HPP
#define PANZER_EVALUATOR_PARAMETER_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {
    
  template <typename EvalT> class ScalarParameterEntry;

//! Constant parameter from sacado parameter library
PHX_EVALUATOR_CLASS(Parameter)
  
  PHX::MDField<ScalarT, Cell, Point> target_field;

  std::size_t cell_data_size;

  ScalarT initial_value;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > param;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_ParameterT.hpp"

#endif
