#ifndef __Panzer_ResponseAggregator_IPCoordinates_T_hpp__
#define __Panzer_ResponseAggregator_IPCoordinates_T_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer {

// useful for cloning and the factory mechanism
template <typename TraitsT>
ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
ResponseAggregator_IPCoordinates() :
  first(true)
{}

template <typename TraitsT>
ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
ResponseAggregator_IPCoordinates(const Teuchos::ParameterList & p) :
  first(true)
{}

template <typename TraitsT>
Teuchos::RCP<ResponseAggregatorBase<TraitsT> > 
ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
clone(const Teuchos::ParameterList & p) const
{ 
   return Teuchos::rcp(new ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>(p)); 
}

//! Build response data for a specified set of fields (ResponseAggregator decides data layout)
template <typename TraitsT>
Teuchos::RCP<ResponseData<TraitsT> > 
ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
buildResponseData(const std::vector<std::string> & fields) const
{
   // simply build data object fully allocated and initialized
   Teuchos::RCP<ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> > data  
         = Teuchos::rcp(new ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT>());
   data->allocateAndInitializeData(fields);
   return data;
}

//! Build an evaluator for the set of fields to be aggregated (calculated) together
template <typename TraitsT>
void ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<ResponseData<TraitsT> > & data,
                             const Teuchos::ParameterList & p) const
{
   typedef ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT> ThisType;

   Teuchos::RCP<ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> > ipc_data =
     Teuchos::rcp_dynamic_cast<ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> >(data);

   typename Teuchos::RCP<IPCoordinates<panzer::Traits::Residual,TraitsT> > eval = 
     Teuchos::rcp(new panzer::IPCoordinates<panzer::Traits::Residual,TraitsT>(1,ipc_data->getNonconstCoords()));

   fm.template registerEvaluator<panzer::Traits::Residual>(eval);

   fm.template requireField<panzer::Traits::Residual>(eval->template getEvaluatedField().fieldTag());
}

//! perform global reduction on this set of response data
template <typename TraitsT>
void
ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const
{ }

template <typename TraitsT>
void ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT>::
aggregateResponses(Response<TraitsT> & dest,const std::list<Teuchos::RCP<const Response<TraitsT> > > & sources) const
{
  Teuchos::RCP<std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > block_agg_coords =
    Teuchos::rcp(new std::vector<typename IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate>);

  for (typename std::list<Teuchos::RCP<const Response<TraitsT> > >::const_iterator block = sources.begin();  block != sources.end(); ++block) {
    
    TEUCHOS_ASSERT((*block)->hasParameterList());
    
    Teuchos::RCP<Teuchos::ParameterList> pl = (*block)->getParameterList();

    Teuchos::RCP<std::vector<panzer::IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > coords = pl->get<Teuchos::RCP<std::vector<panzer::IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate> > >("IP Coordinates");

    for (typename std::vector<panzer::IPCoordinates<panzer::Traits::Residual, Traits>::Coordinate>::const_iterator point = coords->begin(); point != coords->end(); ++point) {
      block_agg_coords->push_back(*point);
    }
    
  }

  Teuchos::RCP<Teuchos::ParameterList> block_agg_pl = Teuchos::parameterList();
  block_agg_pl->set("IP Coordinates",block_agg_coords);
  dest.setParameterList(block_agg_pl);

}

}

#endif
