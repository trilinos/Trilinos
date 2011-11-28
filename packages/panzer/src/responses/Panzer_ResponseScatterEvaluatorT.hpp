#ifndef __Panzer_ResponseScatterEvaluatorT_hpp__
#define __Panzer_ResponseScatterEvaluatorT_hpp__

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <boost/io/ios_state.hpp>
#include <iomanip>

namespace panzer {

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
ResponseScatterEvaluator<EvalT,Traits,AggregatorT>
::ResponseScatterEvaluator(const std::string & name,
                           const Teuchos::RCP<panzer::ResponseData<Traits> > & data,
                           const Teuchos::RCP<const AggregatorT> & aggregator,
                           const std::vector<std::string> & responseNames,
                           int worksetSize)
{
  using Teuchos::RCP;

  // read in some information from response parameter list
  responseData_ = data;
  responseAggregator_ = aggregator;

  // build dummy tag to register with field manager
  responseDummyTag_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>("Response Scatter: " + name,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // add evaluated dummy filed and dependent response fields
  this->addEvaluatedField(*responseDummyTag_);

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  for (std::vector<std::string>::const_iterator s=responseNames.begin(); s!=responseNames.end(); ++s) {
    PHX::MDField<ScalarT,Cell> field(*s,dl_cell);
    responseFields_.push_back(field); // record field in evaluator
    this->addDependentField(field);   // regiseter this as a dependent field
  }
   
  // add dependent fields to evaluator

  std::string n = "Response Scatter: " + name;
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
void ResponseScatterEvaluator<EvalT,Traits,AggregatorT>::
postRegistrationSetup(typename Traits::SetupData sd,PHX::FieldManager<Traits>& fm)
{
   for(typename std::vector<PHX::MDField<ScalarT,Cell> >::iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field)
     this->utils.setFieldData(*field,fm);
}

//**********************************************************************
template <typename EvalT,typename Traits,typename AggregatorT>
void ResponseScatterEvaluator<EvalT,Traits,AggregatorT>::
evaluateFields(typename Traits::EvalData workset)
{ 
  if (workset.num_cells == 0)
    return;

  // do some work
  responseAggregator_->template evaluateFields(workset,*responseData_,responseFields_);
}

} // namespace panzer

#endif
