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

#ifndef __Panzer_ResponseFunctional_Aggregator_impl_hpp__
#define __Panzer_ResponseFunctional_Aggregator_impl_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_CommHelpers.hpp"
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

//! Build an evaluator ofr the set of fields to be aggregated (calculated) together
template <typename TraitsT>
void ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<ResponseData<TraitsT> > & data,
                             const PhysicsBlock & pb,
                             const Teuchos::ParameterList & p) const
{
   typedef ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT> ThisType;

   // build useful evaluator
   Teuchos::RCP<PHX::Evaluator<TraitsT> > eval = Teuchos::rcp(
         new ResponseScatterEvaluator<panzer::Traits::Residual,TraitsT,ThisType>("Functional Response",
                                                                        data,
                                                                        Teuchos::rcpFromRef(*this),
                                                                        data->getFields(),
                                                                        p.get<int>("Workset Size")));

   // add and require fields from aggregator constructed evaluator
   fm.template registerEvaluator<panzer::Traits::Residual>(eval);
   for(std::size_t i=0;i<eval->evaluatedFields().size();i++)
      fm.template requireField<panzer::Traits::Residual>(*(eval->evaluatedFields()[i]));
}

//! Aggregate fields into a specific data object
template <typename TraitsT>
template <typename FieldT>
void ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>::
evaluateFields(panzer::Workset & wkst,ResponseData<TraitsT> & in_data,
               const std::vector<FieldT> & fields) const
{
   ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> & data 
         = Teuchos::dyn_cast<ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> >(in_data); // dynamic cast to correct data type

   std::vector<typename TraitsT::RealType> & dataVec = data.getData();

   TEUCHOS_ASSERT(fields.size()==dataVec.size()); // sanity check

   // loop over reponse fields
   for(std::size_t i=0;i<fields.size();i++) {
      const PHX::MDField<panzer::Traits::Residual::ScalarT,Cell> & field = fields[i]; // this also forces type

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

#include "Panzer_ResponseFunctional_AggregatorSG_impl.hpp"

#endif
