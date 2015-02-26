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

#include <boost/unordered_map.hpp>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Panzer_config.hpp"
#include "Panzer_ResponseBase.hpp"
#include "Panzer_FieldManagerBuilder.hpp"

#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_TypeAssocMap.hpp"

#include "Panzer_ResponseEvaluatorFactory_TemplateManager.hpp"

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

   ResponseLibrary();

   /** Build an initialized response library. By default this
     * method does not initialize the response library to be a residual
     * type. This can be set at runtime to build only residual responses
     * by setting the <code>residualType</code> argument to true.
     */
   ResponseLibrary(const Teuchos::RCP<WorksetContainer> & wc,
                   const Teuchos::RCP<const UniqueGlobalIndexerBase> & ugi,
                   const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof,
                   bool residualType=false); 

   ResponseLibrary(const ResponseLibrary & rl);

   /** Initialize the response library with the appropriate objects.
     */
   void initialize(const Teuchos::RCP<WorksetContainer> & wc,
                   const Teuchos::RCP<const UniqueGlobalIndexerBase> & ugi,
                   const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof); 

   /** Initialize the response library with the appropriate objects. This is
     * in the case that no respones will be added an only a residual is 
     * desired. If <code>addResponse</code> is called then an exception will
     * be raised. 
     */
   void initializeResidualType(const Teuchos::RCP<WorksetContainer> & wc,
                               const Teuchos::RCP<const UniqueGlobalIndexerBase> & ugi,
                               const Teuchos::RCP<const LinearObjFactory<TraitsT> > & lof); 


   /** Initialize the response library from a previously construct response library.
     */
   void initialize(const ResponseLibrary & rl);

   /** Get the internally stored workset container, note this is non-const because
     * the workset container is mostly a non-const object (uses lots of lazy evaluation).
     */
   Teuchos::RCP<WorksetContainer> getWorksetContainer() const
   { return wkstContainer_; }

   //! Get the internally stored global indexer
   Teuchos::RCP<const UniqueGlobalIndexerBase> getGlobalIndexer() const 
   { return globalIndexer_; }

   //! Get the internally stored linear object factory
   Teuchos::RCP<const LinearObjFactory<TraitsT> > getLinearObjFactory() const 
   { return linObjFactory_; }

   /** Add a volumetric response using the response factory builder.
     *
     * \param[in] responseName Name of the response to be added.
     * \param[in] blocks Element blocks to evaluate the response over
     * \param[in] builder Builder that builds the correct response object.
     */
   template <typename ResponseEvaluatorFactory_BuilderT>
   void addResponse(const std::string & responseName,
                    const std::vector<std::string> & blocks,
                    const ResponseEvaluatorFactory_BuilderT & builder); 

   /** Add a surface response using the response factory builder.
     *
     * \param[in] responseName Name of the response to be added.
     * \param[in] sideset_blocks Side set and element blocks to evaluate the response over
     *                           (sideset name is first followed by element block id)
     * \param[in] builder Builder that builds the correct response object.
     */
   template <typename ResponseEvaluatorFactory_BuilderT>
   void addResponse(const std::string & responseName,
                    const std::vector<std::pair<std::string,std::string> > & sideset_blocks,
                    const ResponseEvaluatorFactory_BuilderT & builder); 

   /** Add a response specified by a list of WorksetDescriptor objects. The specifics of the
     * response are specified by the response factory builder.
     *
     * \param[in] responseName Name of the response to be added.
     * \param[in] wkst_desc A vector of descriptors describing the types of elements
     *                                that make up the response.
     * \param[in] builder Builder that builds the correct response object.
     */
   template <typename ResponseEvaluatorFactory_BuilderT>
   void addResponse(const std::string & responseName,
                    const std::vector<WorksetDescriptor> & wkst_desc,
                    const ResponseEvaluatorFactory_BuilderT & builder); 
                   
   /** Access a response by name and evaluation type.
     *
     * \param[in] responseName Name of the response to be retrieved.
     *
     * \return Returns a nonnull response object if it exists, otherwise
     *         it returns null.
     */
   template <typename EvalT>
   Teuchos::RCP<ResponseBase> getResponse(const std::string & responseName) const;

   /** Get the set of responses corresponding to a particular evaluation type. This will
     * overwrite (<code>clear</code>) the vector.
     *
     * \param[in,out] responses Vector over the responses, the responses know their own names!
     */
   template <typename EvalT>
   void getResponses(std::vector<Teuchos::RCP<ResponseBase> > & responses) const;

   /** Setup up field managers for all responses. Once this method is called
     * no other responses can be added. An exception is thrown if they are.
     */
   void buildResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file=false,
         const std::string& graphviz_file_prefix="")
   { buildResponseEvaluators(physicsBlocks,Teuchos::null,cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix); }

   /** Setup up field managers for all responses. Once this method is called
     * no other responses can be added. An exception is thrown if they are.
     */
   void buildResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const panzer::EquationSetFactory & eqset_factory,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file=false,
         const std::string& graphviz_file_prefix="")
   { buildResponseEvaluators(physicsBlocks,Teuchos::ptrFromRef(eqset_factory),cm_factory,closure_models,user_data,write_graphviz_file,graphviz_file_prefix); }

   /** Setup up field managers for a residual response. This method can only be called
     * if the residual response has been setup.
     */
   void buildResidualResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const panzer::EquationSetFactory & eqset_factory,
         const std::vector<BC> & bcs,
         const panzer::BCStrategyFactory & bc_factory,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file=false,
         const std::string& graphviz_file_prefix="");

   /** Have the response evaluators been built? True only if 
     * <code>buildResponseEvaluators</code> has been called and run to completion.
     */ 
   bool responseEvaluatorsBuilt() const
   { return responseEvaluatorsBuilt_; }

   /** Add response objects to assembly data. 
     */
   template <typename EvalT> 
   void addResponsesToInArgs(panzer::AssemblyEngineInArgs & input_args) const;

   /** Evaluate response library for a particular evaluation type.
     */
   template <typename EvalT> 
   void evaluate(const panzer::AssemblyEngineInArgs& input_args);

   /** Print the contents of this response library.
     */
   void print(std::ostream & os) const;

   void useClosureModelByEBlockInResponse(bool value)
   { closureModelByEBlock_ = value; }

   void disableGather(bool value)
   { disableGather_ = value; }

   void disableScatter(bool value)
   { disableScatter_ = value; }

   bool isResidualType() const 
   { return residualType_; }

protected:

   /** Setup up field managers for all responses. Once this method is called
     * no other responses can be added. An exception is thrown if they are.
     */
   void buildResponseEvaluators(
         const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
         const Teuchos::Ptr<const panzer::EquationSetFactory> & eqset_factory,
         const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
         const Teuchos::ParameterList& closure_models,
         const Teuchos::ParameterList& user_data,
         const bool write_graphviz_file,
         const std::string& graphviz_file_prefix);

   /** Add a residual response.
     */
   void addResidualResponse();

   //! A struct for handling function overloading
   template <typename EvalT> struct Overloader {};

   /** Add in the residual responses to the input arguments. Note only residual and Jacobian
     * calls currently work!
     */
   void addResidualResponsesToInArgs(Overloader<typename TraitsT::Residual>,panzer::AssemblyEngineInArgs & input_args) const;

   /** Add in the residual responses to the input arguments. Note only residual and Jacobian
     * calls currently work!
     */
   void addResidualResponsesToInArgs(Overloader<typename TraitsT::Jacobian>,panzer::AssemblyEngineInArgs & input_args) const;

private:

   Teuchos::RCP<WorksetContainer> wkstContainer_;
   Teuchos::RCP<const UniqueGlobalIndexerBase> globalIndexer_;
   Teuchos::RCP<const LinearObjFactory<TraitsT> > linObjFactory_;

   typedef TypeAssocMap<panzer::Traits::EvalTypes,Teuchos::RCP<ResponseBase> > Response_TemplateManager;

   Teuchos::RCP<FieldManagerBuilder> fmb2_;
   AssemblyEngine_TemplateManager<panzer::Traits> ae_tm2_;

   typedef boost::unordered_map<panzer::BC,
                                Teuchos::RCP<std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > >,
                                BC::BCHash,BC::BCEquality > BCHashMap;

   // Store up response factories by element block
   boost::unordered_map<WorksetDescriptor,
                        std::vector<std::pair<std::string,Teuchos::RCP<ResponseEvaluatorFactory_TemplateManager<TraitsT> > > > > respFactories_;
   BCHashMap respBCFactories_;
   std::size_t nextBC_id;

   //! Store all the response objects 
   boost::unordered_map<std::string, Response_TemplateManager> responseObjects_;
   bool closureModelByEBlock_;
   bool disableGather_;
   bool disableScatter_;
   bool residualType_;

   bool responseEvaluatorsBuilt_;

   struct Printer {
     const Response_TemplateManager & tm_;
     std::ostream & os_;
     Printer(const Response_TemplateManager & tm,std::ostream & os) : tm_(tm), os_(os) {}
     template <typename T> void operator()(T) const { 
//       os_ << PHX::TypeString<T>::value << "=";
       if(tm_.get<T>()!=Teuchos::null) 
         os_ << "ON ";
       else
         os_ << "OFF ";
     }
  };
};

}

#include "Panzer_ResponseLibrary_impl.hpp"

#endif
