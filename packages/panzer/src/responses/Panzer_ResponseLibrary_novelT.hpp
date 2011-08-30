#ifndef __Panzer_ResponseLibrary_novelT_hpp__
#define __Panzer_ResponseLibrary_novelT_hpp__

#include "Panzer_ResponseContainer_novel.hpp"

namespace panzer {
namespace novel {

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
buildVolumeFieldManagersFromResponses(const Teuchos::ParameterList & pl)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
  
   // for each element block, build and register evaluators for all required field managers
   typename std::map<std::string,Teuchos::RCP<RespContVector> >::iterator itr;
   for(itr=rsvdVolResp_.begin();itr!=rsvdVolResp_.end();++itr) {
      const std::string & eBlock = itr->first;
      RespContVector & contVector = *(itr->second);

      // build anew field manager of same size as contVector (number of EvalTypes)
      RCP<FMVector> fmVector = rcp(new FMVector(contVector.size()));
      for(std::size_t i=0;i<contVector.size();i++) {

         // if not container has been constructed, don't build
         // a field manager
         if(contVector[i]==Teuchos::null) 
            continue;

         // build and register new field manager
         (*fmVector)[i] = rcp(new PHX::FieldManager<TraitsT>);
         contVector[i]->registerResponses(*(*fmVector)[i],pl);
      }

      volFieldManagers_[eBlock] = fmVector;
   }
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
evaluateVolumeFieldManagers()
{
   // fm.preEvaluate<EvalT>(0);
   // fm.evaluateFields<EvalT>(*workset);
   // fm.postEvaluate<EvalT>(0);
}

}
}

#endif
