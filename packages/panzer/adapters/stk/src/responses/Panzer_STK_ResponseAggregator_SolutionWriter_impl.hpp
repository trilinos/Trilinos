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

#ifndef __Panzer_ResponseAggregator_SolutionWriter_impl_hpp__
#define __Panzer_ResponseAggregator_SolutionWriter_impl_hpp__

#include "Panzer_config.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ScatterFields.hpp"
#include "Panzer_STK_ScatterVectorFields.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"

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
   using Teuchos::RCP;
   using Teuchos::rcp;

   typedef ResponseAggregator_SolutionWriter<EvalT,TraitsT> ThisType;

   const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases = pb.getBases();
   std::map<std::string,std::vector<std::string> > basisBucket;
   bucketByBasisType(pb.getProvidedDOFs(),basisBucket);

   // add this for HCURL and HDIV basis, only want to add them once: evaluate vector fields at centroid
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   RCP<panzer::PointRule> centroidRule;
   for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator itr=bases.begin();
       itr!=bases.end();++itr) {

     if(itr->second->isVectorBasis()) {
        centroidRule = rcp(new panzer::PointRule("Centroid",1,pb.cellData()));

        // compute centroid
        Intrepid::FieldContainer<double> centroid;
        computeReferenceCentroid(bases,pb.cellData().baseCellDimension(),centroid);

        // build pointe values evaluator
        RCP<PHX::Evaluator<panzer::Traits> > evaluator  = 
           rcp(new panzer::PointValues_Evaluator<EvalT,TraitsT>(centroidRule,centroid));
        fm.template registerEvaluator<EvalT>(evaluator);

        break; // get out of the loop, only need one evaluator
     }
   }

   // add evaluators for each field
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

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
      else if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
         TEUCHOS_ASSERT(centroidRule!=Teuchos::null);

         // register basis values evaluator
         {
           Teuchos::RCP<PHX::Evaluator<TraitsT> > evaluator  
              = Teuchos::rcp(new panzer::BasisValues_Evaluator<EvalT,TraitsT>(centroidRule,basis));
           fm.template registerEvaluator<EvalT>(evaluator);
         }

         // add a DOF_PointValues for each field
         std::vector<std::string> pointFields;
         for(std::size_t f=0;f<fields.size();f++) {
            Teuchos::ParameterList p;
            p.set("Name",fields[f]);
            p.set("Basis",basis);
            p.set("Point Rule",centroidRule.getConst());
            Teuchos::RCP<PHX::Evaluator<TraitsT> > evaluator  
               = Teuchos::rcp(new panzer::DOF_PointValues<EvalT,TraitsT>(p));

            fm.template registerEvaluator<EvalT>(evaluator);

            pointFields.push_back(fields[f]+"_"+centroidRule->getName());
         }

         // add the scatter field evaluator for this basis
         {
            Teuchos::RCP<PHX::Evaluator<TraitsT> > evaluator  
               = Teuchos::rcp(new panzer_stk::ScatterVectorFields<EvalT,TraitsT>("STK HCURL Scatter Basis " +basis->name(),
                                                                                 mesh_,centroidRule,fields));

            fm.template registerEvaluator<EvalT>(evaluator);
            fm.template requireField<EvalT>(*evaluator->evaluatedFields()[0]); // require the dummy evaluator
         }
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

template <typename EvalT,typename TraitsT>
void ResponseAggregator_SolutionWriter<EvalT,TraitsT>::
computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases,
                         int baseDimension,
                         Intrepid::FieldContainer<double> & centroid) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   centroid.resize(1,baseDimension);

   // loop over each possible basis
   for(std::map<std::string,RCP<panzer::PureBasis> >::const_iterator itr=bases.begin();
       itr!=bases.end();++itr) {

      // see if this basis has coordinates
      RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis = itr->second->getIntrepidBasis();
      RCP<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<double> > > basisCoords 
         = rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<double> > >(intrepidBasis);

      if(basisCoords==Teuchos::null) // no coordinates...move on
         continue;

      // we've got coordinates, lets commpute the "centroid"
      Intrepid::FieldContainer<double> coords(intrepidBasis->getCardinality(),
                                              intrepidBasis->getBaseCellTopology().getDimension());
      basisCoords->getDofCoords(coords);
      TEUCHOS_ASSERT(coords.rank()==2);
      TEUCHOS_ASSERT(coords.dimension(1)==baseDimension);

      for(int i=0;i<coords.dimension(0);i++)
         for(int d=0;d<coords.dimension(1);d++)
            centroid(0,d) += coords(i,d);

      // take the average
      for(int d=0;d<coords.dimension(1);d++)
         centroid(0,d) /= coords.dimension(0);

      return;
   }

   // no centroid was found...die
   TEUCHOS_ASSERT(false);
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
