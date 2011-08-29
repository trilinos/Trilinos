#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_HPP

#include <iostream>
#include <string>

#include "Panzer_config.hpp"
#include "Panzer_Dimension.hpp"

#include "Panzer_ResponseFunctional_Aggregator.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

/** This class handles responses with values aggregated
  * on each finite element cell.
  *
  * The input parameter list as the following parameter fields
     "Name": type std::string
     "Response Data": type RCP<panzer::ResponseData<Traits> >
     "Response Aggregator": type RCP<panzer::ResponseAggregator<Traits> >
     "Response Names": type RCP<const std::vector<std::string> >
     "Workset Size": type int
  */
PHX_EVALUATOR_CLASS(ResponseScatterEvaluator)
  
  Teuchos::RCP<PHX::FieldTag> responseDummyTag_;
  Teuchos::RCP<panzer::ResponseFunctional_Data<EvalT,Traits> > responseData_;
  Teuchos::RCP<const panzer::ResponseFunctional_Aggregator<EvalT,Traits> > responseAggregator_;
  std::vector<PHX::MDField<ScalarT,Cell> > responseFields_;
    
PHX_EVALUATOR_CLASS_END

}

#include "Panzer_ResponseScatterEvaluatorT.hpp"

#endif
