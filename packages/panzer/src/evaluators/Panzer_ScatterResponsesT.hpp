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

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(ScatterResponses,p)
{
  // read in some information from response parameter list
  comm_ = p.get< Teuchos::RCP<const Teuchos::Comm<int> > >("Comm");
  std::vector<std::string> names = *p.get<Teuchos::RCP<std::vector<std::string> > >("Names");
  int worksetSize = p.get<int>("Workset Size");

  // build response fields to look up
  responseFields_.clear();
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  for (typename std::vector<std::string>::const_iterator name = names.begin(); name != names.end(); ++name)
    responseFields_.push_back(PHX::MDField<ScalarT,Cell>(*name, dl_cell));

  // allocate space for the response values
  responseValues_.push_back(responseFields_.size());
  global_responseValues_.push_back(responseFields_.size());

  // build dummy tag to register with field manager
  responseDummyTag_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>("Response Scatter",Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // add evaluated dummy filed and dependent response fields
  this->addEvaluatedField(*responseDummyTag_);
  for (typename std::vector<PHX::MDField<ScalarT,Cell> >::const_iterator field = responseFields_.begin();
       field != responseFields_.end(); ++field) {
    this->addDependentField(*field);
  }

  std::string n = "Responses";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(ScatterResponses,sd,fm)
{
  for(typename std::vector<PHX::MDField<ScalarT,Cell> >::iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field)
    this->utils.setFieldData(*field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(ScatterResponses,workset)
{ 
  if (workset.num_cells == 0)
    return;

  std::size_t fieldIndex=0;
  for(typename std::vector<PHX::MDField<ScalarT,Cell> >::const_iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field, ++fieldIndex) {
    
    // could be replace with combine mode template
    for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
      responseValues_[fieldIndex] += Sacado::ScalarValue<ScalarT>::eval((*field)(cell)); 
  }
}

//**********************************************************************
PHX_PRE_EVALUATE_FIELDS(ScatterResponses,data)
{
}

//**********************************************************************
PHX_POST_EVALUATE_FIELDS(ScatterResponses,data)
{
  // do communication
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, static_cast<int>(responseValues_.size()),
                     &responseValues_[0], &global_responseValues_[0]);
}

} // namespace panzer

#endif

