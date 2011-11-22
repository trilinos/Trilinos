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
template <typename EvalT,typename Traits>
ResponseScatterEvaluator<EvalT,Traits>
::ResponseScatterEvaluator(const std::string & name,
                           const Teuchos::RCP<panzer::ResponseData<Traits> > & data,
                           const Teuchos::RCP<const panzer::ResponseAggregator<EvalT,Traits> > & aggregator,
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

PHX_EVALUATOR_CTOR(ResponseScatterEvaluator,p)
{
   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                              "Please do not call the ResponseScatterEvalautor(Teuchos::ParameterList) constructor. Use "
                              "another constructor instead!");                     
/*
  using Teuchos::RCP;

  // read in some information from response parameter list
  std::string name = p.get<std::string>("Name");
  responseData_ = p.get<RCP<panzer::ResponseData<Traits> > >("Response Data");
  responseAggregator_ = p.get<RCP<const panzer::ResponseAggregator<ScalarT,Traits> > >("Response Aggregator");
  const std::vector<std::string> & names = *p.get<RCP<const std::vector<std::string> > >("Response Names");
  int worksetSize = p.get<int>("Workset Size");

  // build dummy tag to register with field manager
  responseDummyTag_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>("Response Scatter: " + name,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // add evaluated dummy filed and dependent response fields
  this->addEvaluatedField(*responseDummyTag_);

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  for (std::vector<std::string>::const_iterator name = names.begin(); name != names.end(); ++name) {
    PHX::MDField<ScalarT,Cell> field(*name,dl_cell);
    responseFields_.push_back(field); // record field in evaluator
    this->addDependentField(field);   // regiseter this as a dependent field
  }
   
  // add dependent fields to evaluator

  std::string n = "Response Scatter: " + name;
  this->setName(n);
*/
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(ResponseScatterEvaluator,sd,fm)
{
   for(typename std::vector<PHX::MDField<ScalarT,Cell> >::iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field)
     this->utils.setFieldData(*field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(ResponseScatterEvaluator,workset)
{ 
  if (workset.num_cells == 0)
    return;

  // do some work
  responseAggregator_->evaluateFields(workset,*responseData_,responseFields_);
}

} // namespace panzer

#endif
