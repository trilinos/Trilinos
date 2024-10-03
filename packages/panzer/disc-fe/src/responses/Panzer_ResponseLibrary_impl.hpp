// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseLibrary_impl_hpp__
#define __Panzer_ResponseLibrary_impl_hpp__

#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_ResponseFactory_BCStrategyAdapter.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Panzer_Response_Residual.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include <unordered_set>

namespace panzer {

template <typename TraitsT>
ResponseLibrary<TraitsT>::ResponseLibrary()
   : nextBC_id(0), closureModelByEBlock_(false), disableGather_(false)
   , disableScatter_(true), residualType_(false), responseEvaluatorsBuilt_(false) 
{
}

template <typename TraitsT>
ResponseLibrary<TraitsT>::ResponseLibrary(const Teuchos::RCP<WorksetContainer> & wc,
                                          const Teuchos::RCP<const GlobalIndexer> & ugi,
                                          const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof,
                                          bool residualType)
   : nextBC_id(0), closureModelByEBlock_(false), disableGather_(false)
   , disableScatter_(true), residualType_(false), responseEvaluatorsBuilt_(false)
{
  if(residualType)
    initializeResidualType(wc,ugi,lof);
  else 
    initialize(wc,ugi,lof);
}

template <typename TraitsT>
ResponseLibrary<TraitsT>::ResponseLibrary(const ResponseLibrary<TraitsT> & rl)
   : nextBC_id(0), closureModelByEBlock_(false), disableGather_(false)
   , disableScatter_(true), residualType_(false), responseEvaluatorsBuilt_(false)
{
  initialize(rl);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
initialize(const Teuchos::RCP<WorksetContainer> & wc,
           const Teuchos::RCP<const GlobalIndexer> & ugi,
           const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof)
{
   disableScatter_ = true;
   residualType_ = false;

   wkstContainer_ = wc;
   globalIndexer_ = ugi;
   linObjFactory_ = lof;
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
initializeResidualType(const Teuchos::RCP<WorksetContainer> & wc,
                       const Teuchos::RCP<const GlobalIndexer> & ugi,
                       const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof)
{
   disableScatter_ = false; // we want equation set scatters for this
                            // residual type response
   residualType_ = true;

   wkstContainer_ = wc;
   globalIndexer_ = ugi;
   linObjFactory_ = lof;

   // add the response reponse object
   addResidualResponse();
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
initialize(const ResponseLibrary<TraitsT> & rl)
{
   if(rl.residualType_)
     initializeResidualType(rl.wkstContainer_,rl.globalIndexer_,rl.linObjFactory_);
   else
     initialize(rl.wkstContainer_,rl.globalIndexer_,rl.linObjFactory_);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
copyResponses(const ResponseLibrary & rl)
{
  TEUCHOS_ASSERT(false);
}

namespace panzer_tmp {
  // This is a builder for building a ResponseBase object by evaluation type
  template <typename TraitsT>
  class ResponseBase_Builder {
    Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > respFact_;
    std::string respName_;
    std::vector<WorksetDescriptor> wkstDesc_;

  public:

    // ResponseBase_Builder(const Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > & respFact,
    //                      const std::string & respName, const std::vector<std::string> & eBlocks)
    //   : respFact_(respFact), respName_(respName)
    // {
    //   for(std::size_t i=0;i<eBlocks.size();i++)
    //     wkstDesc_.push_back(blockDescriptor(eBlocks[i]));
    // }

    ResponseBase_Builder(const Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > & respFact,
                         const std::string & respName, const std::vector<std::pair<std::string,std::string> > & sidesets)
      : respFact_(respFact), respName_(respName)
    {
      for(std::size_t i=0;i<sidesets.size();i++)
        wkstDesc_.push_back(sidesetDescriptor(sidesets[i].first,sidesets[i].second));
    }

    ResponseBase_Builder(const Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > & respFact,
                         const std::string & respName, const std::vector<WorksetDescriptor> & wkstDesc)
      : respFact_(respFact), respName_(respName), wkstDesc_(wkstDesc)
    { }

    template <typename T>
    Teuchos::RCP<ResponseBase> build() const 
    { 
      Teuchos::RCP<const panzer::ResponseEvaluatorFactoryBase> baseObj = respFact_->template getAsBase<T>();
     
      // only build this templated set of objects if the there is something to build them with
      if(baseObj!=Teuchos::null && baseObj->typeSupported()) {
        Teuchos::RCP<ResponseBase> resp = baseObj->buildResponseObject(respName_,wkstDesc_); 

        return resp;
      }

      return Teuchos::null;
    }
  };

  // This is a builder for building a ResponseBase object by evaluation type
  template <typename TraitsT>
  class ResidualResponse_Builder {
    std::string respName_;
    Teuchos::RCP<const LinearObjFactory<TraitsT> > lof_;

  public:

    ResidualResponse_Builder(const std::string & respName,const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof)
      : respName_(respName), lof_(lof)
    { }

    template <typename T>
    Teuchos::RCP<ResponseBase> build() const 
    { return Teuchos::rcp(new Response_Residual<T>(respName_,lof_)); }
  };
}

template <typename TraitsT>
template <typename ResponseEvaluatorFactory_BuilderT>
void ResponseLibrary<TraitsT>::
addResponse(const std::string & responseName,
            const std::vector<std::string> & blocks,
            const ResponseEvaluatorFactory_BuilderT & builder) 
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   TEUCHOS_TEST_FOR_EXCEPTION(residualType_,std::invalid_argument,
                              "panzer::ResponseLibrary::addResponse: Method can't be called when the "
                              "response library is a \"residualType\"!");

   // build response factory objects for each evaluation type
   RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > modelFact_tm
        = rcp(new ResponseEvaluatorFactory_TemplateManager<TraitsT>);
   modelFact_tm->buildObjects(builder);

   std::vector<WorksetDescriptor> wkst_desc;
   for(std::size_t i=0;i<blocks.size();i++)
      wkst_desc.push_back(blockDescriptor(blocks[i]));

   addResponse(responseName,wkst_desc,modelFact_tm);
}

template <typename TraitsT>
template <typename ResponseEvaluatorFactory_BuilderT>
void ResponseLibrary<TraitsT>::
addResponse(const std::string & responseName,
            const std::vector<std::pair<std::string,std::string> > & sideset_blocks,
            const ResponseEvaluatorFactory_BuilderT & builder) 
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   TEUCHOS_TEST_FOR_EXCEPTION(residualType_,std::invalid_argument,
                              "panzer::ResponseLibrary::addResponse: Method can't be called when the "
                              "response library is a \"residualType\"!");

   // build response factory objects for each evaluation type
   RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > modelFact_tm
        = rcp(new ResponseEvaluatorFactory_TemplateManager<TraitsT>);
   modelFact_tm->buildObjects(builder);

   // build a response object for each evaluation type
   panzer_tmp::ResponseBase_Builder<TraitsT> respData_builder(modelFact_tm,responseName,sideset_blocks);
   responseObjects_[responseName].buildObjects(respData_builder);

   // associate response objects with all element blocks required
   for(std::size_t i=0;i<sideset_blocks.size();i++) {
     std::string sideset = sideset_blocks[i].first;
     std::string blockId = sideset_blocks[i].second;

     BC bc(nextBC_id,BCT_Neumann,sideset,blockId,"Whatever",responseName+"_BCStrategy");

     // allocate the vector for "bc", if it hasn't yet been allocated
     RCP<std::vector<std::pair<std::string,RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > > block_tm
        = respBCFactories_[bc];
     if(block_tm==Teuchos::null) {
       block_tm = Teuchos::rcp(new std::vector<std::pair<std::string,RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > >); 
       respBCFactories_[bc] = block_tm;
     }

     // add response factory TM to vector that stores them
     block_tm->push_back(std::make_pair(responseName,modelFact_tm));
 
     nextBC_id++;
   }
}

template <typename TraitsT>
template <typename ResponseEvaluatorFactory_BuilderT>
void ResponseLibrary<TraitsT>::
addResponse(const std::string & responseName,
            const std::vector<WorksetDescriptor> & wkst_desc,
            const ResponseEvaluatorFactory_BuilderT & builder) 
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION(residualType_,std::invalid_argument,
                             "panzer::ResponseLibrary::addResponse: Method can't be called when the "
                             "response library is a \"residualType\"!");

  if(wkst_desc[0].useSideset() && !wkst_desc[0].sideAssembly()) {
    // this is a simple side integration, use the "other" addResponse method

    std::vector<std::pair<std::string,std::string> > sideset_blocks;
    for(std::size_t i=0;i<wkst_desc.size();i++) {
      std::string sideset = wkst_desc[i].getSideset();
      std::string blockId = wkst_desc[i].getElementBlock();
      sideset_blocks.push_back(std::make_pair(sideset,blockId));
    }

    // add in the response (as a side set)
    addResponse(responseName,sideset_blocks,builder);

    return;
  }

  // build response factory objects for each evaluation type
  RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > modelFact_tm
       = rcp(new ResponseEvaluatorFactory_TemplateManager<TraitsT>);
  modelFact_tm->buildObjects(builder);

  addResponse(responseName,wkst_desc,modelFact_tm);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResponse(const std::string & responseName,
            const std::vector<WorksetDescriptor> & wkst_desc,
            const Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > & modelFact_tm)
{
  // build a response object for each evaluation type
  panzer_tmp::ResponseBase_Builder<TraitsT> respData_builder(modelFact_tm,responseName,wkst_desc);
  responseObjects_[responseName].buildObjects(respData_builder);

  // associate response objects with all workset descriptors
  for(std::size_t i=0;i<wkst_desc.size();i++) {
    const WorksetDescriptor & desc = wkst_desc[i];

    // add response factory TM to vector that stores them
    respFactories_[desc].push_back(std::make_pair(responseName,modelFact_tm));
  }
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResidualResponse()
{
   std::string responseName = "RESIDUAL";

   // setup responses to be constructed
   panzer_tmp::ResidualResponse_Builder<TraitsT> respData_builder(responseName,linObjFactory_);

   // build all the response objects (for each evaluation type)
   responseObjects_[responseName].buildObjects(respData_builder);
}

template <typename TraitsT>
template <typename EvalT>
Teuchos::RCP<ResponseBase> ResponseLibrary<TraitsT>::
getResponse(const std::string & responseName) const
{
   typedef std::unordered_map<std::string, Response_TemplateManager> HashMap;
   HashMap::const_iterator itr = responseObjects_.find(responseName);

   // response was not in list of responses
   if(itr==responseObjects_.end())
     return Teuchos::null;

   // response was found, return it
   return itr->second.get<EvalT>();
}

template <typename TraitsT>
template <typename EvalT>
void ResponseLibrary<TraitsT>::
getResponses(std::vector<Teuchos::RCP<ResponseBase> > & responses) const
{
   typedef std::unordered_map<std::string, Response_TemplateManager> HashMap;

   responses.clear();

   // loop over all respones adding them to the vector
   for(HashMap::const_iterator itr=responseObjects_.begin();itr!=responseObjects_.end();++itr)
     responses.push_back(itr->second.get<EvalT>());
}

template <typename TraitsT>
class RVEF2 : public GenericEvaluatorFactory {
public:
   typedef std::unordered_map<WorksetDescriptor,
                                std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > > RespFactoryTable;

   RVEF2(const Teuchos::ParameterList & userData,RespFactoryTable & rft)
     : userData_(userData), rft_(rft) {}

   bool registerEvaluators(PHX::FieldManager<TraitsT> & fm,const WorksetDescriptor & wd, const PhysicsBlock & pb) const
   {
     using Teuchos::RCP;
     using Teuchos::rcp;

     TEUCHOS_ASSERT(wd.getElementBlock()==pb.elementBlockID());

     // because we reduce the physics blocks to only ones we need, this find should succeed
     typename RespFactoryTable::iterator itr=rft_.find(wd);

     TEUCHOS_ASSERT(itr!=rft_.end() && itr->second.size()>0);

     // loop over each template manager in the block, and then each evaluation type
     std::vector<std::pair<std::string,RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > & respFacts = itr->second;
     for(std::size_t i=0;i<respFacts.size();i++) {
       std::string responseName = respFacts[i].first;
       RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > fact = respFacts[i].second;
       // what is going on here?
       if(fact==Teuchos::null)
         continue;
 
       // loop over each evaluation type
       for(typename ResponseEvaluatorFactory_TemplateManager<TraitsT>::iterator rf_itr=fact->begin();
           rf_itr!=fact->end();++rf_itr) {

         // not setup for this template type, ignore it
         if(rf_itr.rcp()==Teuchos::null || !rf_itr.rcp()->typeSupported()) 
           continue;

         // build and register evaluators, store field tag, make it required
         rf_itr->buildAndRegisterEvaluators(responseName,fm,pb,userData_);     
       }
     }
    
     return true;
   }

private:
  const Teuchos::ParameterList & userData_;
  RespFactoryTable & rft_;
};

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
buildResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const Teuchos::Ptr<const panzer::EquationSetFactory> & eqset_factory,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file,
         const std::string& graphviz_file_prefix)
{
   using Teuchos::RCP;

   TEUCHOS_TEST_FOR_EXCEPTION(residualType_,std::invalid_argument,
                              "panzer::ResponseLibrary::buildResponseEvaluators: Method can't be called when the "
                              "response library is a \"residualType\"!");

   typedef std::unordered_map<WorksetDescriptor,
                                std::vector<std::pair<std::string,RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > > RespFactoryTable;

   // first compute subset of physics blocks required to build responses
   ////////////////////////////////////////////////////////////////////////////////

   std::vector<Teuchos::RCP<panzer::PhysicsBlock> > requiredVolPhysicsBlocks;
   std::vector<WorksetDescriptor> requiredWorksetDesc;
   for(typename RespFactoryTable::const_iterator itr=respFactories_.begin();
       itr!=respFactories_.end();++itr) {
     // is there something to do?
     if(itr->second.size()==0) 
       continue;

     const WorksetDescriptor & wd = itr->first;

     // find physics block with right element block
     bool failure = true;
     for(std::size_t i=0;i<physicsBlocks.size();i++) {
       if(physicsBlocks[i]->elementBlockID()==wd.getElementBlock()) {
         requiredVolPhysicsBlocks.push_back(physicsBlocks[i]);
         failure = false;
         break;
       } 
     }

     // we must find at least one physics block
     // TEUCHOS_ASSERT(!failure);
     if(!failure)
       requiredWorksetDesc.push_back(wd);
   }

   // build boundary response array
   std::vector<BC> bcs;
   for(typename BCHashMap::const_iterator itr=respBCFactories_.begin();
       itr!=respBCFactories_.end();++itr)
     bcs.push_back(itr->first);

   // second construct generic evaluator factory from required response factories
   ////////////////////////////////////////////////////////////////////////////////
   RVEF2<TraitsT> rvef2(user_data,respFactories_);

   // third build field manager builder using the required physics blocks
   ////////////////////////////////////////////////////////////////////////////////

   response_bc_adapters::BCFactoryResponse bc_factory(respBCFactories_);

   // don't build scatter evaluators
   fmb2_ = Teuchos::rcp(new FieldManagerBuilder(disableScatter_,disableGather_)); 

   fmb2_->setWorksetContainer(wkstContainer_);
   fmb2_->setupVolumeFieldManagers(requiredVolPhysicsBlocks,requiredWorksetDesc,cm_factory,closure_models,*linObjFactory_,user_data,rvef2,closureModelByEBlock_);
   if(eqset_factory==Teuchos::null)
     fmb2_->setupBCFieldManagers(bcs,physicsBlocks,cm_factory,bc_factory,closure_models,*linObjFactory_,user_data);
   else
     fmb2_->setupBCFieldManagers(bcs,physicsBlocks,*eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory_,user_data);

   if(write_graphviz_file) {
     fmb2_->writeVolumeGraphvizDependencyFiles("Response_Volume_"+graphviz_file_prefix,requiredVolPhysicsBlocks);
     fmb2_->writeBCGraphvizDependencyFiles("Response_Surface_"+graphviz_file_prefix);
   }

   // fourth build assembly engine from FMB
   ////////////////////////////////////////////////////////////////////////////////

   AssemblyEngine_TemplateBuilder builder(fmb2_,linObjFactory_); 
   ae_tm2_.buildObjects(builder);

   responseEvaluatorsBuilt_ = true;
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
buildResidualResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const panzer::EquationSetFactory & eqset_factory,
         const std::vector<BC> & bcs,
         const panzer::BCStrategyFactory& bc_factory,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file,
         const std::string& graphviz_file_prefix)
{
   using Teuchos::RCP;

   TEUCHOS_TEST_FOR_EXCEPTION(!residualType_,std::invalid_argument,
                              "panzer::ResponseLibrary::buildResidualResponseEvaluators: Method can only be called when the "
                              "response library is a \"residualType\"!");

   // don't build scatter evaluators
   fmb2_ = Teuchos::rcp(new FieldManagerBuilder(disableScatter_,disableGather_)); 

   fmb2_->setWorksetContainer(wkstContainer_);
   fmb2_->setupVolumeFieldManagers(physicsBlocks,cm_factory,closure_models,*linObjFactory_,user_data);
   fmb2_->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,cm_factory,bc_factory,closure_models,*linObjFactory_,user_data);

   // Print Phalanx DAGs
   if (write_graphviz_file){
     fmb2_->writeVolumeGraphvizDependencyFiles("ResidualResponse_Volume_"+graphviz_file_prefix,physicsBlocks);
     fmb2_->writeBCGraphvizDependencyFiles("ResidualResponse_Surface_"+graphviz_file_prefix);
   }

   // fourth build assembly engine from FMB
   ////////////////////////////////////////////////////////////////////////////////

   AssemblyEngine_TemplateBuilder builder(fmb2_,linObjFactory_); 
   ae_tm2_.buildObjects(builder);

   responseEvaluatorsBuilt_ = true;
}

template <typename TraitsT>
template <typename EvalT> 
void ResponseLibrary<TraitsT>::
addResponsesToInArgs(panzer::AssemblyEngineInArgs & input_args) const
{
  std::vector<Teuchos::RCP<ResponseBase> > responses;
  this->getResponses<EvalT>(responses);

  // add all responses to input args  
  if(!residualType_) {
    for(std::size_t i=0;i<responses.size();i++) {
      if(responses[i]!=Teuchos::null) {
        input_args.addGlobalEvaluationData(responses[i]->getLookupName(),responses[i]);
      }
    }
  }
  else { // residualType_ == true
    addResidualResponsesToInArgs(Overloader<EvalT>(),input_args);
  }
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResidualResponsesToInArgs(Overloader<typename TraitsT::Residual>,panzer::AssemblyEngineInArgs & input_args) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef typename TraitsT::Residual EvalT;
  typedef typename TraitsT::RealType ScalarT;

  // extract the residual response
  RCP<Response_Residual<EvalT> > resp = rcp_dynamic_cast<Response_Residual<EvalT> >(getResponse<EvalT>("RESIDUAL"));
  resp->initializeResponse();

  // setup the local ghosted container
  if(ghostedContainer_==Teuchos::null)
    ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

  // replace ghosted container with local one
  const RCP<ThyraObjContainer<ScalarT> > thGhostedContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(ghostedContainer_);
  input_args.ghostedContainer_ = ghostedContainer_;
  
  // convert responses into thyra object
  const RCP<ThyraObjContainer<ScalarT> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(input_args.container_);

  // set the ghosted and unique residual
  thGhostedContainer->set_f_th(resp->getGhostedResidual());
  thGlobalContainer->set_f_th(resp->getResidual());

  TEUCHOS_ASSERT(thGhostedContainer->get_f_th()!=Teuchos::null);
  TEUCHOS_ASSERT(thGlobalContainer->get_f_th()!=Teuchos::null);

  // clear out ghosted residual
  Thyra::assign(thGhostedContainer->get_f_th().ptr(),0.0);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResidualResponsesToInArgs(Overloader<typename TraitsT::Jacobian>,panzer::AssemblyEngineInArgs & input_args) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef typename TraitsT::Jacobian EvalT;
  typedef typename TraitsT::RealType ScalarT;

  // extract the residual response
  RCP<Response_Residual<EvalT> > resp = rcp_dynamic_cast<Response_Residual<EvalT> >(getResponse<EvalT>("RESIDUAL"));
  resp->initializeResponse();

  // setup the local ghosted container
  if(ghostedContainer_==Teuchos::null)
    ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

  // replace ghosted container with local one
  const RCP<ThyraObjContainer<ScalarT> > thGhostedContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(ghostedContainer_);
  input_args.ghostedContainer_ = ghostedContainer_;

  // convert responses into thyra object
  const RCP<ThyraObjContainer<ScalarT> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(input_args.container_);

  // set the ghosted and unique residual
  thGhostedContainer->set_A_th(resp->getGhostedJacobian());

  RCP<Thyra::VectorBase<ScalarT> > dummy_f = Thyra::createMember(resp->getJacobian()->range());
  thGlobalContainer->set_f_th(dummy_f);
  thGlobalContainer->set_A_th(resp->getJacobian());

  // Zero values in ghosted container objects
  thGhostedContainer->initializeMatrix(0.0);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResidualResponsesToInArgs(Overloader<typename TraitsT::Tangent>,panzer::AssemblyEngineInArgs & input_args) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef typename TraitsT::Tangent EvalT;
  typedef typename TraitsT::RealType ScalarT;

  // extract the residual response
  RCP<Response_Residual<EvalT> > resp = rcp_dynamic_cast<Response_Residual<EvalT> >(getResponse<EvalT>("RESIDUAL"));
  resp->initializeResponse();

  // setup the local ghosted container
  if(ghostedContainer_==Teuchos::null)
    ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

  // replace ghosted container with local one
  const RCP<ThyraObjContainer<ScalarT> > thGhostedContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(ghostedContainer_);
  input_args.ghostedContainer_ = ghostedContainer_;

  // convert responses into thyra object
  const RCP<ThyraObjContainer<ScalarT> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(input_args.container_);

  // At this point it isn't clear what to do for Tangent.  We probably need to extend the linear object containers
  // to support df/dp.

  /*
  // set the ghosted and unique residual
  thGhostedContainer->set_A_th(resp->getGhostedJacobian());

  RCP<Thyra::VectorBase<ScalarT> > dummy_f = Thyra::createMember(resp->getJacobian()->range());
  thGlobalContainer->set_f_th(dummy_f);
  thGlobalContainer->set_A_th(resp->getJacobian());

  // Zero values in ghosted container objects
  thGhostedContainer->initializeMatrix(0.0);
  */
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template <typename TraitsT>
void ResponseLibrary<TraitsT>::
addResidualResponsesToInArgs(Overloader<typename TraitsT::Hessian>,panzer::AssemblyEngineInArgs & input_args) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef typename TraitsT::Hessian EvalT;
  typedef typename TraitsT::RealType ScalarT;

  // extract the residual response
  RCP<Response_Residual<EvalT> > resp = rcp_dynamic_cast<Response_Residual<EvalT> >(getResponse<EvalT>("RESIDUAL"));
  resp->initializeResponse();

  // setup the local ghosted container
  if(ghostedContainer_==Teuchos::null)
    ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

  // replace ghosted container with local one
  const RCP<ThyraObjContainer<ScalarT> > thGhostedContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(ghostedContainer_);
  input_args.ghostedContainer_ = ghostedContainer_;

  // convert responses into thyra object
  const RCP<ThyraObjContainer<ScalarT> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<ThyraObjContainer<ScalarT> >(input_args.container_);

  // set the ghosted and unique residual
  thGhostedContainer->set_A_th(resp->getGhostedHessian());

  RCP<Thyra::VectorBase<ScalarT> > dummy_f = Thyra::createMember(resp->getHessian()->range());
  thGlobalContainer->set_f_th(dummy_f);
  thGlobalContainer->set_A_th(resp->getHessian());

  // Zero values in ghosted container objects
  thGhostedContainer->initializeMatrix(0.0);
}
#endif

template <typename TraitsT>
template <typename EvalT> 
void ResponseLibrary<TraitsT>::
evaluate(const panzer::AssemblyEngineInArgs& input_args)
{
   ae_tm2_.template getAsObject<EvalT>()->evaluate(input_args);
}

template <typename TraitsT>
void ResponseLibrary<TraitsT>::
print(std::ostream & os) const
{
   typedef std::unordered_map<std::string, Response_TemplateManager> RespObjType;

   for(RespObjType::const_iterator itr=responseObjects_.begin();itr!=responseObjects_.end();++itr) { 
     std::string respName = itr->first;
     os << "Response \"" << respName << "\": ";
     Sacado::mpl::for_each<typename Response_TemplateManager::types_vector>(Printer(itr->second,os));
     os << std::endl;
   }
}

}

#endif
