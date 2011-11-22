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
     "Response Aggregator": type RCP<panzer::ResponseAggregator<ScalarT,Traits> >
     "Response Names": type RCP<const std::vector<std::string> >
     "Workset Size": type int
  */
PHX_EVALUATOR_CLASS(ResponseScatterEvaluator)
  
  Teuchos::RCP<PHX::FieldTag> responseDummyTag_;
  Teuchos::RCP<panzer::ResponseData<Traits> > responseData_;
  Teuchos::RCP<const panzer::ResponseAggregator<EvalT,Traits> > responseAggregator_;
  std::vector<PHX::MDField<ScalarT,Cell> > responseFields_;

public:
  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator(const std::string & name,
                           const Teuchos::RCP<panzer::ResponseData<Traits> > & data,
                           const Teuchos::RCP<const panzer::ResponseAggregator<EvalT,Traits> > & aggregator,
                           const std::vector<std::string> & responseNames,
                           int worksetSize);
    
PHX_EVALUATOR_CLASS_END

}

#include "Panzer_ResponseScatterEvaluatorT.hpp"

#endif
