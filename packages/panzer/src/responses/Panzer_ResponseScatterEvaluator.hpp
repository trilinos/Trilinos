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
  */
template<typename EvalT, typename Traits,typename AggregatorT>
class ResponseScatterEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                                 public PHX::EvaluatorDerived<EvalT, Traits>  { 
public:

  //! A constructor with concrete arguments instead of a parameter list.
  ResponseScatterEvaluator(const std::string & name,
                           const Teuchos::RCP<panzer::ResponseData<Traits> > & data,
                           const Teuchos::RCP<const AggregatorT> & aggregator,
                           const std::vector<std::string> & responseNames,
                           int worksetSize);

  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:

  typedef typename EvalT::ScalarT ScalarT;

  Teuchos::RCP<PHX::FieldTag> responseDummyTag_;
  Teuchos::RCP<panzer::ResponseData<Traits> > responseData_;
  Teuchos::RCP<const AggregatorT> responseAggregator_;
  std::vector<PHX::MDField<ScalarT,Cell> > responseFields_;

};

}

#include "Panzer_ResponseScatterEvaluator_impl.hpp"

#endif
