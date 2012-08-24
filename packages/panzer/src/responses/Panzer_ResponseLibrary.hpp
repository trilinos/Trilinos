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

#ifndef __Panzer_ResponseLibrary_hpp__
#define __Panzer_ResponseLibrary_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Response.hpp"
#include "Panzer_ResponseAggregator_Manager.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseFunctional_Aggregator.hpp"
#include "Panzer_RLDynamicDispatch.hpp"

#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

namespace panzer {

template <typename TraitsT> class ResponseContainerBase;
template <typename EvalT,typename TraitsT> class ResponseContainer;

/** This contains, collects and serves as a resource for
  * responses computed by panzer. This functions as a library
  * where there are many "responses" maintained (as many as
  * a user adds).  When a response is maintained that simply
  * means there is a mechansim to "reserve" it. A response is
  * not required by any field manager until a user "reserves"
  * it. The reservation process is done by response name and 
  * the element block or BC it is associated with. The field
  * tag specified in the <code>addResponse</code> is used only
  * for the "name" field, however that use of <code>PHX::FieldTag</code>
  * reminds the user that something better be in the evaluation tree.
  */
template <typename TraitsT>
class ResponseLibrary {
public:
   typedef typename TraitsT::EvalTypes TypeSeq;

   ResponseLibrary() 
   {
      // build dynamic dispatch objects
      dynamicDispatch_.buildObjects(Teuchos::ptrFromRef(*this)); 
   }

   ResponseLibrary(const Teuchos::RCP<WorksetContainer> & wc,
                   const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
                   const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof); 

   ResponseLibrary(const ResponseLibrary & rl);

   /** Initialize the response library with the appropriate objects.
     */
   void initialize(const Teuchos::RCP<WorksetContainer> & wc,
                   const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
                   const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof); 

   /** Initialize the response library from a previously construct response library.
     */
   void initialize(const ResponseLibrary & rl);

   /** Asks, does this string correspond to a response type
     * in this library?
     */
   template <typename EvalT>
   bool isResponseType(const std::string & type) const
   { 
     return getAggregatorManager().template isAggregator<EvalT>(type);
   }

   //! Get a particular aggregator 
   template <typename EvalT>
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type) const
   {
      return getAggregatorManager().template getAggregator<EvalT>(type);
   }

   //! dynamic dispatch version
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type,const std::string & evalType) const
   { return dynamicDispatch_.getAggregator(type,evalType); }

   /** \defgroup volume
     * Volume methods for volumetric reponses
     * @{
     */

   /** User access to a volume response object. This will throw if called before
     * the volume field managers (all blocks) are evaluated.
     */
   Teuchos::RCP<const Response<TraitsT> > getVolumeResponse(const ResponseId & rid,
                                                            const std::string & eBlock) const;

   /** Get a particular volume response by label.
     */ 
   Teuchos::RCP<const Response<TraitsT> > getBlockAggregatedVolumeResponseByLabel(const std::string & label) const;

   //! Reserve a response for actual calculation (by response id and element block).
   template <typename EvalT>
   void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock);

   /** Reserve a response for calculation (by response id, element block, and evalution type).
     * The evaluation type must be a string that is define by <code>PHX::TypeString<EvalT>==evalType</code>
     * a member of the <code>Traits::EvalTypes</code> type list.
     */
   void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType);

   /** Reserve a labeled response over a set of element blocks and evaluation types. A labeled
     * response is a first class citizen, and the expected mechanism used to access and register
     * responses. 
     *
     * \param[in] label User readable label for the response
     * \param[in] rid Response identifier containing the relevant field and type of response.
     * \param[in] eBlocks Element blocks to be aggregated over.
     * \param[in] evalTypes String of evaluation types for the response (Residual,Jacobian, etc...)
     */
   void reserveLabeledBlockAggregatedVolumeResponse(const std::string & label,const ResponseId & rid,
						    const std::list<std::string> & eBlocks,
						    const std::list<std::string> & evalTypes);

   /** Veryify that this response and element block are actual valid choices
     * for the evaluation type. This is optional error checking but makes debugging
     * simplier.
     */
   template <typename EvalT>
   bool validateResponseIdInElementBlock(const ResponseId & rid,const std::string & eBlock) const 
   { return true; }

   /** This method builds the volume field managers from the reserved
     * responses. It also registers a number of evaluators, using the closure
     * model and equation set factories. Unlike in the assembly engine only gather
     * evaluators, DOF, and Gradient evaluators are automatically included. This
     * method also stores the worksets passed in to be used when the evaluate method
     * is called.
     */
   void buildVolumeFieldManagersFromResponses(
                        const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                        const Teuchos::ParameterList& closure_models,
                        const Teuchos::ParameterList& user_data,
                        const bool write_graphviz_file=false,
                        const std::string& graphviz_file_prefix="");

   /** Evaluate all the volume field managers of a particular evaluator type.
     */
   template <typename EvalT>
   void evaluateVolumeFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
                                    const panzer::AssemblyEngineInArgs & ae_in,
                                    const Teuchos::Comm<int> & comm);

   /** Evaluate all the volume field managers of a particular evaluator type.
     */
   template <typename EvalT>
   void evaluateVolumeFieldManagers(const panzer::AssemblyEngineInArgs & ae_in,
                                    const Teuchos::Comm<int> & comm);

   /** @} */

   /** Returns the set of element blocks required by the ResponseLibrary to 
     * build any responses. This only details those element blocks the library is
     * currently aware of, if new responses are reserved in new element blocks
     * those will not be included. Additionally, this defines the set of worksets
     * required for the evaluating the responses.
     */
   void getRequiredElementBlocks(std::vector<std::string> & eBlocks) const;

   //! Write out all volume containers to a stream
   void printVolumeContainers(std::ostream & os) const;

   //! Access the response aggregator manager
   ResponseAggregator_Manager<TraitsT> & getAggregatorManager()
   { return respAggManager_; }

   //! Access the response aggregator manager
   const ResponseAggregator_Manager<TraitsT> & getAggregatorManager() const
   { return respAggManager_; }

   //! Define some default aggregators used in panzer
   void defineDefaultAggregators()
   { ResponseAggregator_Manager<TraitsT>::defineDefaultAggregators(getAggregatorManager()); }

   //! How many responses are contained by this library.
   std::size_t getLabeledResponseCount() const 
   { return labeledResponses_.size(); }

   //! get all labeled respones
   void getLabeledVolumeResponses(std::vector<Teuchos::RCP<const Response<TraitsT> > > & responses) const;

   //! get volume response labels
   void getVolumeResponseLabels(std::vector<std::string> & labels) const;

   //! Reinitialize the reponse data
   void reinitializeResponseData();

protected:
   //! Access a container field for a specified element block
   template <typename EvalT>
   Teuchos::RCP<ResponseContainerBase<TraitsT> > getVolumeContainer(const std::string & eBlock);

private:
   // This could be a template manager, but this turns out to be more in line with
   // what is desired here. Direct access to the vector.
   typedef std::vector<Teuchos::RCP<ResponseContainerBase<TraitsT> > > RespContVector;
   typedef std::vector<Teuchos::RCP<PHX::FieldManager<TraitsT> > > FMVector;

   struct ResponseDescriptor {
      ResponseId rid;
      std::list<std::string> elmtBlocks;
      std::list<std::string> evalTypes;
   };

   ResponseAggregator_Manager<TraitsT> respAggManager_;

   std::map<std::string,Teuchos::RCP<RespContVector> > rsvdVolResp_;
   std::map<std::string,Teuchos::RCP<PHX::FieldManager<TraitsT> > > volFieldManagers_;

   RLDynamicDispatch<TraitsT> dynamicDispatch_;

   Teuchos::RCP<WorksetContainer> wkstContainer_;

   //! Maps from a "label"->(ResponseId,List of Element Blocks)
   std::map<std::string,ResponseDescriptor> labeledResponses_;

   Teuchos::RCP<UniqueGlobalIndexerBase> globalIndexer_;
   Teuchos::RCP<LinearObjFactory<TraitsT> > linObjFactory_;
};

}

#include "Panzer_ResponseLibrary_impl.hpp"

#endif
