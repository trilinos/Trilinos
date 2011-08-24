#include "Panzer_ResponseContainer.hpp"

#include "Panzer_config.hpp"
#include "Panzer_ScatterResponses.hpp"

namespace panzer {

// Value type responses
//////////////////////////////////////////

ResponseContainer<panzer::Traits::Value>::ResponseContainer()
   : responseFieldsBuilt_(false)
{
}

void ResponseContainer<panzer::Traits::Value>::
registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                  PHX::FieldManager<panzer::Traits> & fm) 
{
   typedef RespType::EvalType EvalType;

   using Teuchos::RCP;
   using Teuchos::rcp;

   comm_ = comm;

   // also build response fields
   if(!responseFieldsBuilt_)
      buildResponseFields(worksetSize);

   std::vector<std::string> names;
   this->getReserved(names);

   // add scatter response
   {
      Teuchos::ParameterList p;
      p.set("Name", getResponseType());
      p.set<Teuchos::RCP<ResponseContainerBase> >("Response Container", Teuchos::rcpFromRef(*this));
 
      RCP<panzer::ScatterResponses<EvalType,panzer::Traits> > e 
         = rcp(new panzer::ScatterResponses<EvalType,panzer::Traits>(p));

      fm.requireField<EvalType>(e->getRequiredFieldTag());
      fm.registerEvaluator<EvalType>(e);
   }

   responseVector_.resize(names.size());
   localResponseVector_.resize(names.size());
}

void ResponseContainer<panzer::Traits::Value>::
addDependentFields(PHX::EvaluatorWithBaseImpl<panzer::Traits> & eval) const
{
  for(std::vector<PHX::MDField<RespType::ScalarT,Cell> >::const_iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field) 
    eval.addDependentField(*field);
}

void ResponseContainer<panzer::Traits::Value>::
buildResponseFields(int worksetSize)
{
  std::vector<std::string> names;
  this->getReserved(names);

  // build response fields to look up
  responseFields_.clear();
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  for (std::vector<std::string>::const_iterator name = names.begin(); name != names.end(); ++name)
    responseFields_.push_back(PHX::MDField<RespType::ScalarT,Cell>(*name, dl_cell));

  responseFieldsBuilt_ = true;
}

void ResponseContainer<panzer::Traits::Value>::
setFieldData(PHX::FieldManager<panzer::Traits> & fm)
{
   localResponseVector_.clear(); // can sum local responses
   localResponseVector_.resize(responseVector_.size());
 
   // load up
   std::vector<PHX::MDField<RespType::ScalarT,Cell> >::iterator itr;
   for(itr=responseFields_.begin();itr!=responseFields_.end();++itr) 
      fm.getFieldData<RespType::ScalarT,RespType::EvalType>(*itr);
}

void ResponseContainer<panzer::Traits::Value>::
aggregateFieldsLocally(const panzer::Workset & wkst)
{
  std::size_t fieldIndex=0;
  for(std::vector<PHX::MDField<RespType::ScalarT,Cell> >::const_iterator field = responseFields_.begin();
      field != responseFields_.end(); ++field, ++fieldIndex) {
    
    // value does a summation operation
    for (std::size_t cell = 0; cell < wkst.num_cells; ++cell)
      localResponseVector_[fieldIndex] += Sacado::ScalarValue<RespType::ScalarT>::eval((*field)(cell)); 
  }
}

void ResponseContainer<panzer::Traits::Value>::
reduceFieldsGlobally()
{
  // do communication
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, static_cast<int>(responseVector_.size()),
                     &localResponseVector_[0], &responseVector_[0]);
}

//! This returns the response for a given field
Teuchos::RCP<const Response<panzer::Traits::Value> > ResponseContainer<panzer::Traits::Value>::
getResponse(const std::string & name) const
{
   std::size_t i = this->getFieldIndex(name); 
   return Teuchos::rcp(new Response<RespType>(responseVector_[i])); 
}

// Derivative type responses
//////////////////////////////////////////

ResponseContainer<panzer::Traits::Derivative>::ResponseContainer()
{
}

}
