#ifndef PANZER_SCATTER_RESPONSES_HPP
#define PANZER_SCATTER_RESPONSES_HPP

#include <iostream>
#include <string>
#include "Panzer_config.hpp"
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Comm.hpp"

#include "Sacado_Traits.hpp"

namespace panzer {
    
/** This class handles responses with values aggregated
  * on each finite element cell.
  */
PHX_EVALUATOR_CLASS_PP(ScatterResponses)
  
  std::vector<PHX::MDField<ScalarT,Cell> > responseFields_;
  Teuchos::RCP<PHX::FieldTag> responseDummyTag_;
    
  std::vector<typename Sacado::ScalarType<ScalarT>::type > responseValues_;
  std::vector<typename Sacado::ScalarType<ScalarT>::type > global_responseValues_;

  Teuchos::RCP<const Teuchos::Comm<int> > comm_;

public:
  const PHX::FieldTag& getRequiredFieldTag()
  { return *responseDummyTag_; }


PHX_EVALUATOR_CLASS_END

}

#include "Panzer_ScatterResponsesT.hpp"

#endif
