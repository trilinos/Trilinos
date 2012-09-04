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

#ifndef __Panzer_ResponseLibrary_impl_hpp__
#define __Panzer_ResponseLibrary_impl_hpp__

#include "Panzer_ResponseContainer.hpp"

namespace panzer {

template <typename TraitsT>
ResponseLibrary<TraitsT>::ResponseLibrary(const Teuchos::RCP<WorksetContainer> & wc,
                                          const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
                                          const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof)
   : respAggManager_(ugi,lof), wkstContainer_(wc), globalIndexer_(ugi), linObjFactory_(lof)
{
   // build dynamic dispatch objects
   dynamicDispatch_.buildObjects(Teuchos::ptrFromRef(*this)); 
}

template <typename TraitsT>
ResponseLibrary<TraitsT>::ResponseLibrary(const ResponseLibrary<TraitsT> & rl)
   : respAggManager_(rl.globalIndexer_,rl.linObjFactory_), wkstContainer_(rl.wkstContainer_)
   , globalIndexer_(rl.globalIndexer_), linObjFactory_(rl.linObjFactory_)
{
   // build dynamic dispatch objects
   dynamicDispatch_.buildObjects(Teuchos::ptrFromRef(*this)); 
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
initialize(const Teuchos::RCP<WorksetContainer> & wc,
           const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
           const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof)
{
   respAggManager_ .initialize(ugi,lof);
   wkstContainer_ = wc;
   globalIndexer_ = ugi;
   linObjFactory_ = lof;
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
initialize(const ResponseLibrary<TraitsT> & rl)
{
   respAggManager_.initialize(rl.globalIndexer_,rl.linObjFactory_);
   wkstContainer_ = rl.wkstContainer_;
   globalIndexer_ = rl.globalIndexer_; 
   linObjFactory_ = rl.linObjFactory_;
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock)
{
   validateResponseIdInElementBlock<EvalT>(rid,eBlock);

   int idx = Sacado::mpl::find<TypeSeq,EvalT>::value;
   int sz = Sacado::mpl::size<TypeSeq>::value;

   // response container vector for this element block does not yet
   // exist, build and initialize it!
   Teuchos::RCP<RespContVector> respContMngr = rsvdVolResp_[eBlock];
   if(respContMngr==Teuchos::null) {
      respContMngr = Teuchos::rcp(new RespContVector(sz,Teuchos::null));
      rsvdVolResp_[eBlock] = respContMngr;
   }

   // if container does not yet exist, build and initialize it
   Teuchos::RCP<ResponseContainerBase<TraitsT> > container = (*respContMngr)[idx];
   if(container==Teuchos::null) {
      container = Teuchos::rcp(new ResponseContainer<EvalT,TraitsT>);
      container->setResponseLibrary(Teuchos::rcpFromRef(*this));
      (*respContMngr)[idx] = container;
   } 

   // reserve this respoinse id
   container->reserve(rid);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType)
{
   dynamicDispatch_.reserveVolumeResponse(rid,eBlock,evalType);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
reserveLabeledBlockAggregatedVolumeResponse(const std::string & label,const ResponseId & rid,
					    const std::list<std::string> & eBlocks,
					    const std::list<std::string> & evalTypes)
{
   TEUCHOS_TEST_FOR_EXCEPTION(labeledResponses_.find(label)!=labeledResponses_.end(),std::logic_error,
                      "ResponseLibrary::reserveLabeledVolumeResponse: Adding response labeled \""+label+"\" "
                      "failed because response label has already been added!");

   // add labeled responses
   labeledResponses_[label].rid = rid;
   labeledResponses_[label].elmtBlocks = eBlocks;
   labeledResponses_[label].evalTypes = evalTypes;

   // loop over element blocks
   for(std::list<std::string>::const_iterator eBlk=eBlocks.begin(); 
       eBlk!=eBlocks.end();++eBlk) {
      // loop over evaluation types
      for(std::list<std::string>::const_iterator eType=evalTypes.begin(); 
          eType!=evalTypes.end();++eType) {
 
         // reserve this response
         reserveVolumeResponse(rid,*eBlk,*eType);
      }
   }
}

template <typename TraitsT>
template <typename EvalT>
Teuchos::RCP<ResponseContainerBase<TraitsT> > ResponseLibrary<TraitsT>::
getVolumeContainer(const std::string & eBlock)
{
   // validateElementBlock(eBlock); // should I add this?

   int idx = Sacado::mpl::find<TypeSeq,EvalT>::value;
   return rsvdVolResp_[eBlock][idx];
}

template <typename TraitsT>
Teuchos::RCP<const Response<TraitsT> > ResponseLibrary<TraitsT>::
getVolumeResponse(const ResponseId & rid,const std::string & eBlock) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // grab vector of containers associated with the element blocks
   RCP<RespContVector> vec;
   { 
      typename std::map<std::string,Teuchos::RCP<RespContVector> >::const_iterator 
            itr = rsvdVolResp_.find(eBlock);
      TEUCHOS_TEST_FOR_EXCEPTION(itr==rsvdVolResp_.end(),std::logic_error, 
                         "Could not find element block \""+eBlock+"\" in response library");
      vec = itr->second;
   }
   
   // Loop over all reponse containers extracting data and aggregating them into the
   // response
   bool responseDataFound = false;
   RCP<Response<TraitsT> > response = rcp(new Response<TraitsT>(rid)); 
   for(typename RespContVector::const_iterator itr=vec->begin();itr!=vec->end();++itr) {
      // if no container is associated with this evaluation type then move on.
      if(*itr==Teuchos::null) 
         continue; 

      RCP<ResponseData<TraitsT> > data = (*itr)->getResponseData(rid.type); 

      if(data!=Teuchos::null) { // is there an aggregator ("type") for this evaluation type?
         data->fillResponse(rid.name,*response);
         responseDataFound = true;
      }
   }

   TEUCHOS_TEST_FOR_EXCEPTION(!responseDataFound,std::logic_error,
                      "ReponseLibrary::getVolumeResponse could not find any such response \""+rid.getString() +"\""
                      " in element block \""+eBlock+"\"");

   return response;
}

/** Get a particular volume response by label.
  */ 
template <typename TraitsT>
Teuchos::RCP<const Response<TraitsT> > ResponseLibrary<TraitsT>::
getBlockAggregatedVolumeResponseByLabel(const std::string & label) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   typename std::map<std::string,ResponseDescriptor>::const_iterator itr=labeledResponses_.find(label);
   TEUCHOS_TEST_FOR_EXCEPTION(itr==labeledResponses_.end(),std::logic_error,
                      "ResponseLibrary::getVolumeResponseByLabel: Cannot find response labeled \""+label+"\"!");
   
   const ResponseId & rid = itr->second.rid;
   const std::list<std::string> & eBlocks = itr->second.elmtBlocks;
   const std::list<std::string> & evalTypes = itr->second.evalTypes;

   // get responses for each element block
   std::list<RCP<const Response<TraitsT> > > blkResponses;
   for(std::list<std::string>::const_iterator eblkItr=eBlocks.begin();
       eblkItr!=eBlocks.end();++eblkItr) 
      blkResponses.push_back(getVolumeResponse(rid,*eblkItr));

   TEUCHOS_TEST_FOR_EXCEPTION(blkResponses.size()==0,std::logic_error,
                      "ReponseLibrary::getVolumeResponseByLabel: Could not find any response in "
                      "subcontainers for Response label \""+label+"\"!");

   // for each evaluation type use an aggregator to aggregate responses
   RCP<Response<TraitsT> > response = rcp(new Response<TraitsT>(rid)); 
   for(std::list<std::string>::const_iterator eTypeItr=evalTypes.begin();
       eTypeItr!=evalTypes.end();++eTypeItr) {
      getAggregator(rid.type,*eTypeItr).aggregateResponses(*response,blkResponses);
   }

   return response;
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
getRequiredElementBlocks(std::vector<std::string> & eBlocks) const
{
   eBlocks.clear();
   for(typename std::map<std::string,Teuchos::RCP<RespContVector> >::const_iterator itr=rsvdVolResp_.begin();
       itr!=rsvdVolResp_.end();++itr) {
      eBlocks.push_back(itr->first);
   }
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
buildVolumeFieldManagersFromResponses(
                        const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                        const Teuchos::ParameterList& closure_models,
                        const Teuchos::ParameterList& user_data,
                        const bool write_graphviz_file,
                        const std::string& graphviz_file_prefix)
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
   for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
      RCP<panzer::PhysicsBlock> pb = *blkItr;
      std::string blockId = pb->elementBlockID();

      // verify that block is relevant
      typename std::map<std::string,Teuchos::RCP<RespContVector> >::iterator contItr = rsvdVolResp_.find(blockId);
      if(contItr==rsvdVolResp_.end())
         continue;
      RespContVector & contVector = *contItr->second;
  
      // build a field manager object
      Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
            = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
      
      // Choose model sublist for this element block
      std::string response_model_name = "Response Model";
     
      Teuchos::ParameterList tmp_user_data(user_data);
      tmp_user_data.set<bool>("Ignore Scatter",true);
  
      // use the physics block to register evaluators
      pb->buildAndRegisterEquationSetEvaluators(*fm, user_data);
      pb->buildAndRegisterGatherAndOrientationEvaluators(*fm,*linObjFactory_, tmp_user_data);
      pb->buildAndRegisterDOFProjectionsToIPEvaluators(*fm, tmp_user_data);
      pb->buildAndRegisterScatterEvaluators(*fm,*linObjFactory_, tmp_user_data);
      pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, closure_models, user_data);
      if(closure_models.isSublist(response_model_name))
         pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, response_model_name, closure_models, user_data);

      for(std::size_t i=0;i<contVector.size();i++) {

         // if container has not been constructed, don'ta register responses
         if(contVector[i]==Teuchos::null) 
            continue;

         // build and register new field manager
         contVector[i]->registerResponses(*fm,*pb,user_data);
      }
  
      // build the setup data using passed in information
      Traits::SetupData setupData;
      setupData.worksets_ = wkstContainer_->getVolumeWorksets(blockId);
  
      fm->postRegistrationSetup(setupData);
      
      if (write_graphviz_file)
        fm->writeGraphvizFile(graphviz_file_prefix+"Response_"+blockId);

      volFieldManagers_[blockId] = fm;
   }
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
evaluateVolumeFieldManagers(const panzer::AssemblyEngineInArgs & ae_in,
                            const Teuchos::Comm<int> & comm)
{
   int idx = Sacado::mpl::find<TypeSeq,EvalT>::value;

   GlobalEvaluationDataContainer preEvalData;
   preEvalData.addDataObject("Solution Gather Container",ae_in.ghostedContainer_);
   ae_in.fillGlobalEvaluationDataContainer(preEvalData);

   typedef panzer::LinearObjContainer LOC;
   linObjFactory_->globalToGhostContainer(*(ae_in.container_),*(ae_in.ghostedContainer_),LOC::X | LOC::DxDt);

   // std::map<std::string,Teuchos::RCP<PHX::FieldManager<TraitsT> > >::iterator fm_itr;
   typename std::map<std::string,Teuchos::RCP<PHX::FieldManager<TraitsT> > >::iterator fm_itr;
   for(fm_itr=volFieldManagers_.begin();fm_itr!=volFieldManagers_.end();fm_itr++) {
     const std::string & eBlock = fm_itr->first;
     Teuchos::RCP< PHX::FieldManager<TraitsT> > fm = fm_itr->second;
 
     fm->template preEvaluate<EvalT>(preEvalData);

     // loop over all worksets
     for(std::vector<panzer::Workset>::iterator wkst_itr=wkstContainer_->begin(eBlock);
         wkst_itr!=wkstContainer_->end(eBlock);++wkst_itr) {
 
        panzer::Workset& workset = *wkst_itr;
    
        workset.ghostedLinContainer = ae_in.ghostedContainer_;
        workset.linContainer = ae_in.container_;
        workset.alpha = ae_in.alpha;
        workset.beta = ae_in.beta;
        workset.time = ae_in.time;
        workset.evaluate_transient_terms = ae_in.evaluate_transient_terms;
 
        fm->template evaluateFields<EvalT>(workset);
     }
 
     fm->template postEvaluate<EvalT>(0);
 
     // perform global communication
     const RespContVector & contVector = *rsvdVolResp_.find(eBlock)->second;
    
     // if not container has been constructed, don't build
     // a field manager
     if(contVector[idx]!=Teuchos::null) 
        contVector[idx]->globalReduction(comm);
   }
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
evaluateVolumeFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
                            const panzer::AssemblyEngineInArgs & ae_in,
                            const Teuchos::Comm<int> & comm)
{
  int idx = Sacado::mpl::find<TypeSeq,EvalT>::value;

  GlobalEvaluationDataContainer preEvalData;
  preEvalData.addDataObject("Solution Gather Container",ae_in.ghostedContainer_);
  ae_in.fillGlobalEvaluationDataContainer(preEvalData);

  typedef panzer::LinearObjContainer LOC;
  linObjFactory_->globalToGhostContainer(*(ae_in.container_),*(ae_in.ghostedContainer_),LOC::X | LOC::DxDt);

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >::const_iterator itr;
  for(itr=worksets.begin();itr!=worksets.end();++itr) {
    const std::string & eBlock = itr->first;
    std::vector<panzer::Workset> & w = *(itr->second);

    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = volFieldManagers_[eBlock];

    if(fm!=Teuchos::null) {
       fm->preEvaluate<EvalT>(preEvalData);
   
       // Loop over worksets in this element block
       for (std::size_t i = 0; i < w.size(); ++i) {
         panzer::Workset& workset = w[i];
   
         workset.ghostedLinContainer = ae_in.ghostedContainer_;
         workset.linContainer = ae_in.container_;
         workset.alpha = ae_in.alpha;
         workset.beta = ae_in.beta;
         workset.time = ae_in.time;
         workset.evaluate_transient_terms = ae_in.evaluate_transient_terms;

         fm->evaluateFields<EvalT>(workset);
       }
   
       fm->postEvaluate<EvalT>(0);
   
       // perform global communication
       const RespContVector & contVector = *rsvdVolResp_.find(eBlock)->second;
   
       // if not container has been constructed, don't build
       // a field manager
       if(contVector[idx]!=Teuchos::null) 
          contVector[idx]->globalReduction(comm);
    }
  }

}

//! Write out all volume containers to a stream
template <typename TraitsT>
void ResponseLibrary<TraitsT>::
printVolumeContainers(std::ostream & os) const
{
   // loop over all active containers
   for(typename std::map<std::string,Teuchos::RCP<RespContVector> >::const_iterator itr=rsvdVolResp_.begin();
       itr!=rsvdVolResp_.end();++itr) {
      const std::string & eBlock = itr->first;
      const RespContVector & respContVec = *itr->second;

      os << "Element Block = \"" << eBlock << "\"" << std::endl;
      for(std::size_t i=0;i<respContVec.size();i++) {
         if(respContVec[i]!=Teuchos::null) 
            os << "   " << *respContVec[i] << std::endl;
      }
   }
}

//! get all labeled respones
template <typename TraitsT>
void ResponseLibrary<TraitsT>::
getLabeledVolumeResponses(std::vector<Teuchos::RCP<const Response<TraitsT> > > & responses) const
{
   responses.clear();

   for(typename std::map<std::string,ResponseDescriptor>::const_iterator itr=labeledResponses_.begin();
       itr!=labeledResponses_.end();++itr)
      responses.push_back(getBlockAggregatedVolumeResponseByLabel(itr->first));
}

//! get all labeled respones
template <typename TraitsT>
void ResponseLibrary<TraitsT>::
getVolumeResponseLabels(std::vector<std::string> & labels) const
{
   labels.clear();

   for(typename std::map<std::string,ResponseDescriptor>::const_iterator itr=labeledResponses_.begin();
       itr!=labeledResponses_.end();++itr)
      labels.push_back(itr->first);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
reinitializeResponseData() 
{
   // loop over all active containers
   for(typename std::map<std::string,Teuchos::RCP<RespContVector> >::iterator itr=rsvdVolResp_.begin();
       itr!=rsvdVolResp_.end();++itr) {
      RespContVector & respContVec = *itr->second;
      for(std::size_t i=0;i<respContVec.size();i++) {
         if(respContVec[i]!=Teuchos::null)
            respContVec[i]->clear();
      }
   }
}

}

#endif
