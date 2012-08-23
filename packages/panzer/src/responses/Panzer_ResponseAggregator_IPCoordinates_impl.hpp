// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_ResponseAggregator_IPCoordinates_impl_hpp__
#define __Panzer_ResponseAggregator_IPCoordinates_impl_hpp__

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
                             const PhysicsBlock & pb,
                             const Teuchos::ParameterList & p) const
{
   typedef ResponseAggregator_IPCoordinates<panzer::Traits::Residual,TraitsT> ThisType;

   int ip_order = p.sublist("IP Coordinates").get<int>("Integration Order");

   Teuchos::RCP<ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> > ipc_data =
     Teuchos::rcp_dynamic_cast<ResponseAggregator_IPCoordinates_Data<panzer::Traits::Residual,TraitsT> >(data);

   typename Teuchos::RCP<IPCoordinates<panzer::Traits::Residual,TraitsT> > eval = 
     Teuchos::rcp(new panzer::IPCoordinates<panzer::Traits::Residual,TraitsT>(ip_order,ipc_data->getNonconstBlockID(),ipc_data->getNonconstCoords()));

   fm.template registerEvaluator<panzer::Traits::Residual>(eval);

   fm.template requireField<panzer::Traits::Residual>(eval->getEvaluatedField().fieldTag());
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
  Teuchos::RCP< std::map<std::string,Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > > > block_agg_coords =
    Teuchos::rcp(new std::map<std::string,Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > >);

  for (typename std::list<Teuchos::RCP<const Response<TraitsT> > >::const_iterator block = sources.begin();  block != sources.end(); ++block) {
    
    TEUCHOS_ASSERT((*block)->hasParameterList());
    
    Teuchos::RCP<Teuchos::ParameterList> pl = (*block)->getParameterList();

    std::string block_id = pl->get<std::string>("Element Block ID");

    Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > coords = 
      pl->get<Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > >("IP Coordinates");

    (*block_agg_coords)[block_id] = coords;
    
  }

  Teuchos::RCP<Teuchos::ParameterList> block_agg_pl = Teuchos::parameterList();
  block_agg_pl->set("IP Coordinates",block_agg_coords);
  dest.setParameterList(block_agg_pl);

}

}

#endif
