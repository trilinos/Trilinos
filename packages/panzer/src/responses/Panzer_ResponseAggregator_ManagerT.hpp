#ifndef __Panzer_ResponseAggregator_ManagerT_hpp__
#define __Panzer_ResponseAggregator_ManagerT_hpp__

#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

template <typename TraitsT>
template <typename EvalT>
const ResponseAggregatorBase<TraitsT> & ResponseAggregator_Manager<TraitsT>::
getAggregator(const std::string & type) const
{
   Teuchos::RCP<const ResponseAggregatorBase<TraitsT> > agg = getAggregatorManager(type).template getAsBase<EvalT>();

   TEST_FOR_EXCEPTION(agg==Teuchos::null,std::logic_error,
                      "ResponseAggregator_Manager cannot find type \""+type+"\" associated with "
                      "evaluation type \""+PHX::TypeString<EvalT>::value+"\"");

   // guranteed not to fail because of above test
   return *agg;
}

template <typename TraitsT>
template <typename Builder>
void ResponseAggregator_Manager<TraitsT>::
defineAggregatorTypeFromBuilder(const std::string & type,const Builder & builder)
{
   TEST_FOR_EXCEPTION(isAggregator(type),std::logic_error,
                      "ResponeAggregator_Manager: Aggregator of type \""+type+"\" "
                      "already exists! Choose a new type name!");

   Teuchos::RCP<AggregatorManager> aggMngr = Teuchos::rcp(new AggregatorManager);
   aggMngr->buildObjects(builder);
   aggregators_[type] = aggMngr;
} 

//! Is there an aggretator of a particular type.
template <typename TraitsT>
bool ResponseAggregator_Manager<TraitsT>::
isAggregator(const std::string & type) const
{
   return aggregators_.find(type)!=aggregators_.end();
}

 //! Is there an aggregator that uses a particular evaluation type
template <typename TraitsT>
template <typename EvalT>
bool ResponseAggregator_Manager<TraitsT>::
isAggregator() const
{
   // loop over each aggregator type
   for(typename std::map<std::string,Teuchos::RCP<AggregatorManager> >::const_iterator itr = aggregators_.begin();
       itr!=aggregators_.end();++itr) {

      // check if this contains the EvalT evaluation type
      if(isAggregator<EvalT>(itr->first))
         return true;
   }

   return false;
}

 //! Is there an aggregator that is a particular type and uses a specified evaluation type
template <typename TraitsT>
template <typename EvalT>
bool ResponseAggregator_Manager<TraitsT>::
isAggregator(const std::string & type) const
{
   typename std::map<std::string,Teuchos::RCP<AggregatorManager> >::const_iterator itr = aggregators_.find(type);

   // could not find a paricular type
   if(itr==aggregators_.end())
      return false;

   // can we find the evaulation type?
   return (*itr->second).template getAsBase<EvalT>()!=Teuchos::null;
}

template <typename TraitsT>
const typename ResponseAggregator_Manager<TraitsT>::AggregatorManager &
ResponseAggregator_Manager<TraitsT>::
getAggregatorManager(const std::string & type) const
{
   typename std::map<std::string,Teuchos::RCP<AggregatorManager> >::const_iterator itr = aggregators_.find(type);
  
   TEST_FOR_EXCEPTION(itr==aggregators_.end(),std::logic_error,
                      "ResponseAggregator_Manager does not have any aggregators of type \"" + type + "\"!");

   return *(itr->second);
}

template <typename TraitsT>
void ResponseAggregator_Manager<TraitsT>::
defineDefaultAggregators(ResponseAggregator_Manager<TraitsT> & aggMngr)
{
   // add functional aggregator
   {
      ResponseFunctional_Aggregator_Builder builder;
      aggMngr.defineAggregatorTypeFromBuilder("Functional",builder);
   }
}

} // end namespace panzer

#endif
