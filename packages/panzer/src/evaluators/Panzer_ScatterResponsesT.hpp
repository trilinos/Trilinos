#ifndef PANZER_SCATTER_RESPONSES_T_HPP
#define PANZER_SCATTER_RESPONSES_T_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <boost/io/ios_state.hpp>
#include <iomanip>

#include "Panzer_ResponseContainer.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(ScatterResponses,p)
{
  using Teuchos::RCP;

  // read in some information from response parameter list
  std::string name = p.get<std::string>("Name");
  responseContainer_ = p.get<RCP<panzer::ResponseContainerBase> >("Response Container");

  // build dummy tag to register with field manager
  responseDummyTag_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>("Response Scatter: " + name,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // add evaluated dummy filed and dependent response fields
  this->addEvaluatedField(*responseDummyTag_);
   
  // add dependent fields to evaluator
  responseContainer_->addDependentFields(*this);

  std::string n = "Responses: " + name;
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(ScatterResponses,sd,fm)
{
  responseContainer_->setFieldData(fm);

  // for(typename std::vector<PHX::MDField<ScalarT,Cell> >::iterator field = responseFields_.begin();
  //     field != responseFields_.end(); ++field)
  //   this->utils.setFieldData(*field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(ScatterResponses,workset)
{ 
  if (workset.num_cells == 0)
    return;

  responseContainer_->aggregateFieldsLocally(workset);
}

//**********************************************************************
PHX_PRE_EVALUATE_FIELDS(ScatterResponses,data)
{
   responseContainer_->clear();
}

//**********************************************************************
PHX_POST_EVALUATE_FIELDS(ScatterResponses,data)
{
   responseContainer_->reduceFieldsGlobally();
}

} // namespace panzer

#endif
