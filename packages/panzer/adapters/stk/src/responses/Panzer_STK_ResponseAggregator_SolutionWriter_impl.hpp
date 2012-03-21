// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef __Panzer_ResponseAggregator_SolutionWriter_impl_hpp__
#define __Panzer_ResponseAggregator_SolutionWriter_impl_hpp__

#include "Panzer_config.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ScatterFields.hpp"

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer_stk {

// useful for cloning and the factory mechanism
template <typename EvalT,typename TraitsT>
ResponseAggregator_SolutionWriter<EvalT,TraitsT>::
ResponseAggregator_SolutionWriter(const Teuchos::RCP<STK_Interface> & mesh)
   : mesh_(mesh)
{
}

//! Build an evaluator for the set of fields to be aggregated (calculated) together
template <typename EvalT,typename TraitsT>
void ResponseAggregator_SolutionWriter<EvalT,TraitsT>::
registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,const Teuchos::RCP<panzer::ResponseData<TraitsT> > & data,
                             const panzer::PhysicsBlock & pb,
                             const Teuchos::ParameterList & p) const
{
   typedef ResponseAggregator_SolutionWriter<EvalT,TraitsT> ThisType;

   const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases = pb.getBases();
   std::map<std::string,std::vector<std::string> > basisBucket;
   bucketByBasisType(pb.getProvidedDOFs(),basisBucket);

   for(std::map<std::string,std::vector<std::string> >::const_iterator itr=basisBucket.begin();
       itr!=basisBucket.end();++itr) {

      std::string basisName = itr->first;
      const std::vector<std::string> & fields = itr->second;

      std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator found = bases.find(basisName);
      TEUCHOS_ASSERT(found!=bases.end());
      Teuchos::RCP<const panzer::PureBasis> basis = found->second;
      
      // write out nodal fields
      if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {

         Teuchos::RCP<PHX::Evaluator<TraitsT> > eval = 
           Teuchos::rcp(new ScatterFields<EvalT,TraitsT>("STK HGRAD Scatter Basis " +basis->name(),
                                                         mesh_, basis, fields));
   
         // register and require evaluator fields
         fm.template registerEvaluator<EvalT>(eval);
         fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
      }
   }
}

template <typename EvalT,typename TraitsT>
void ResponseAggregator_SolutionWriter<EvalT,TraitsT>::
bucketByBasisType(const std::vector<panzer::StrPureBasisPair> & providedDofs,
                  std::map<std::string,std::vector<std::string> > & basisBucket)
{
   // this should be self explanatory
   for(std::size_t i=0;i<providedDofs.size();i++) {
      std::string fieldName = providedDofs[i].first;
      Teuchos::RCP<const panzer::PureBasis> basis = providedDofs[i].second;

      basisBucket[basis->name()].push_back(fieldName);
   }
}

// Specializations for residual evaluation type

template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > 
ResponseAggregator_SolutionWriter_Builder::
build<panzer::Traits::Residual>() const
{ 
   Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseAggregator_SolutionWriter<panzer::Traits::Residual,panzer::Traits>(mesh_)); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}

#ifdef HAVE_STOKHOS
template < >
Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > 
ResponseAggregator_SolutionWriter_Builder::
build<panzer::Traits::SGResidual>() const
{ 
   Teuchos::RCP<panzer::ResponseAggregatorBase<panzer::Traits> > respAgg =
      Teuchos::rcp(new ResponseAggregator_SolutionWriter<panzer::Traits::SGResidual,panzer::Traits>(mesh_)); 
   respAgg->setLinearObjFactory(getLinearObjFactory());
   respAgg->setGlobalIndexer(getGlobalIndexer());
   return respAgg;
}
#endif

}

#endif
