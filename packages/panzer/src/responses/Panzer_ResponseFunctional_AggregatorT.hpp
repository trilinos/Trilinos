#ifndef __Panzer_ResponseFunctional_AggregatorT_hpp__
#define __Panzer_ResponseFunctional_AggregatorT_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_ResponseScatterEvaluator.hpp"

namespace panzer {

// useful for cloning and the factory mechanism
template <typename TraitsT>
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
ResponseFunctional_Aggregator() 
{}

template <typename TraitsT>
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
ResponseFunctional_Aggregator(const Teuchos::ParameterList & p) 
{}

template <typename TraitsT>
Teuchos::RCP<ResponseAggregatorBase<TraitsT> > 
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
clone(const Teuchos::ParameterList & p) const
{ 
   return Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>(p)); 
}

//! Build response data for a specified set of fields (ResponseAggregator decides data layout)
template <typename TraitsT>
Teuchos::RCP<ResponseData<TraitsT> > 
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
buildResponseData(const std::vector<std::string> & fields) const
{
   // simply build data object fully allocated and initialized
   Teuchos::RCP<ResponseData<TraitsT> > data  
         = Teuchos::rcp(new ResponseFunctional_Data<panzer::Traits::Residual,TraitsT>());
   data->allocateAndInitializeData(fields);
   
   return data;
}

//! Build an evaluator of the set of fields to be aggregated (calculated) together
template <typename TraitsT>
Teuchos::RCP<PHX::Evaluator<TraitsT> > 
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
buildEvaluator(Teuchos::RCP<ResponseData<TraitsT> > & data, const Teuchos::ParameterList & fix_p) const
{
   Teuchos::RCP<ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> > func_data 
      = Teuchos::rcp_dynamic_cast<ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> >(data);

   Teuchos::ParameterList p(fix_p);
   p.set("Name","Functional Response");
   p.set("Response Data",func_data);
   p.set("Response Aggregator",Teuchos::rcpFromRef(*this));
   p.set("Response Names",Teuchos::rcpFromRef(func_data->getFields()));

   p.get<int>("Workset Size");  // check to make sure its there!

   // build useful evaluator
   return Teuchos::rcp(new ResponseScatterEvaluator<panzer::Traits::Residual,panzer::Traits>(p));
}

//! Build an evaluator ofr the set of fields to be aggregated (calculated) together
template <typename TraitsT>
void
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,Teuchos::RCP<ResponseData<TraitsT> > & data,
                             const Teuchos::ParameterList & p) const
{
   Teuchos::RCP<PHX::Evaluator<TraitsT> > eval = buildEvaluator(data,p);

   // add and require fields from aggregator constructed evaluator
   fm.template registerEvaluator<panzer::Traits::Residual>(eval);
   for(std::size_t i=0;i<eval->evaluatedFields().size();i++)
      fm.template requireField<panzer::Traits::Residual>(*(eval->evaluatedFields()[i]));
}

//! Aggregate fields into a specific data object
template <typename TraitsT>
void
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
evaluateFields(panzer::Workset & wkst,ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> & data,
                                           const std::vector<PHX::MDField<panzer::Traits::Residual::ScalarT,Cell> > & fields) const
{
   std::vector<typename TraitsT::RealType> & dataVec = data.getData();

   TEUCHOS_ASSERT(fields.size()==dataVec.size()); // sanity check

   // loop over reponse fields
   for(std::size_t i=0;i<fields.size();i++) {
      const PHX::MDField<panzer::Traits::Residual::ScalarT,Cell> & field = fields[i];

      // loop over cells
      for(std::size_t c=0;c<wkst.num_cells;c++)
         dataVec[i] += field(c);
   }
}

//! perform global reduction on this set of response data
template <typename TraitsT>
void
ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const
{
   std::vector<typename TraitsT::RealType> & dataVec = 
      Teuchos::dyn_cast<ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> >(rd).getData();
   std::vector<typename TraitsT::RealType> dataVec2 = dataVec;

  // do communication
  Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, static_cast<int>(dataVec.size()),
                     &dataVec2[0], &dataVec[0]);
}

template <typename TraitsT>
void ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
aggregateResponses(Response<TraitsT> & dest,const std::list<Teuchos::RCP<const Response<TraitsT> > > & sources) const
{
   typename TraitsT::RealType val = dest.getValue();

   // sum over all values
   typename std::list<Teuchos::RCP<const Response<TraitsT> > >::const_iterator itr;
   for(itr=sources.begin();itr!=sources.end();++itr)
      val += (*itr)->getValue();      

   dest.setValue(val);
}

}

#endif
