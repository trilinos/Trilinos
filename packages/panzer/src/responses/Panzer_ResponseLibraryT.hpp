#ifndef __Panzer_ResponseLibraryT_hpp__
#define __Panzer_ResponseLibraryT_hpp__

#include "Panzer_ResponseContainer.hpp"

namespace panzer {

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
      TEST_FOR_EXCEPTION(itr==rsvdVolResp_.end(),std::logic_error, 
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

   TEST_FOR_EXCEPTION(!responseDataFound,std::logic_error,
                      "ReponseLibrary::getVolumeResponse could not find any such response \""+rid.getString() +"\""
                      " in element block \""+eBlock+"\"");

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
                        const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets,
                        const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                        const Teuchos::ParameterList& ic_block_closure_models,
                        const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
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
      RespContVector & contVector = *rsvdVolResp_.find(blockId)->second;
  
      // build a field manager object
      Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
            = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
      
      // Choose model sublist for this element block
      std::string closure_model_name = "";
      if (ic_block_closure_models.isSublist(blockId))
        closure_model_name = blockId;
      else if (ic_block_closure_models.isSublist("Default"))
        closure_model_name = "Default";
      else 
        TEST_FOR_EXCEPTION(true, std::logic_error, 
                           "Failed to find initial condition for element block \"" << blockId << 
                           "\".  You must provide an initial condition for each element block or set a default!" << ic_block_closure_models);
     
      Teuchos::ParameterList tmp_user_data(user_data);
      tmp_user_data.set<bool>("Ignore Scatter",true);
  
      // use the physics block to register evaluators
      pb->buildAndRegisterEquationSetEvaluators(*fm, user_data);
      pb->buildAndRegisterGatherScatterEvaluators(*fm,lo_factory, tmp_user_data);
      pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, ic_block_closure_models.sublist(closure_model_name), user_data);

      for(std::size_t i=0;i<contVector.size();i++) {

         // if not container has been constructed, don't build
         // a field manager
         if(contVector[i]==Teuchos::null) 
            continue;

         // build and register new field manager
         contVector[i]->registerResponses(*fm,user_data);
      }
  
      // build the setup data using passed in information
      Traits::SetupData setupData;
      setupData.worksets_ = volume_worksets.find(blockId)->second;
  
      fm->postRegistrationSetup(setupData);
      
      if (write_graphviz_file)
        fm->writeGraphvizFile(graphviz_file_prefix+"Response_"+blockId);

      volFieldManagers_[blockId] = fm;
   }
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
evaluateVolumeFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
                            const Teuchos::RCP<panzer::LinearObjContainer> & loc,const Teuchos::Comm<int> & comm)
{
  int idx = Sacado::mpl::find<TypeSeq,EvalT>::value;

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >::const_iterator itr;
  for(itr=worksets.begin();itr!=worksets.end();++itr) {
    const std::string & eBlock = itr->first;
    std::vector<panzer::Workset> & w = *(itr->second);

    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = volFieldManagers_[eBlock];

    if(fm!=Teuchos::null) {
       fm->preEvaluate<EvalT>(0);
   
       // Loop over worksets in this element block
       for (std::size_t i = 0; i < w.size(); ++i) {
         panzer::Workset& workset = w[i];
   
         workset.linContainer = loc;
   
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

}

#endif
