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

#ifndef __Panzer_ResponseAggregator_Manager_hpp__
#define __Panzer_ResponseAggregator_Manager_hpp__

#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include "Phalanx_TemplateManager.hpp"

namespace panzer {

template <typename TraitsT>
class ResponseAggregator_Manager {
public:
   typedef PHX::TemplateManager<typename TraitsT::EvalTypes,ResponseAggregatorBase<TraitsT>,ResponseAggregator<_,TraitsT> >
           AggregatorManager;

   ResponseAggregator_Manager()
      : globalIndexer_(Teuchos::null), linObjFactory_(Teuchos::null) {}

   ResponseAggregator_Manager(const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
                              const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof)
      : globalIndexer_(ugi), linObjFactory_(lof) {}

   void initialize(const Teuchos::RCP<UniqueGlobalIndexerBase> & ugi,
                   const Teuchos::RCP<LinearObjFactory<TraitsT> > & lof)
   { globalIndexer_ = ugi; linObjFactory_ = lof; }

   /** Statically access an aggregator of a particular type
     *
     * \note This method throws no aggregator exists with the 
     *       type string and evaluation type template parameter.
     */
   template <typename EvalT>
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type) const;

   /** Add a user defined aggregator to this container with a particular type.
     * The user must provide a builder as defined in the <code>PHX::TemplateManager</code>
     * documentation. This method will excercise the builder and construct a list of
     * aggregators for each evaluation type. In this case it is safe, and expected,
     * that the builder return Teuchos::null if no aggregator is appropriate for that evalution
     * type. This is not the case for the <code>PHX::TemplateManager</code>
     *
     * \note This method throws if this type has already been defined
     */
   template <typename Builder>
   void defineAggregatorTypeFromBuilder(const std::string & type,const Builder & builder);

   /** \defgroup isAggregator
     * Look ups for if an aggregator is contained in this manager.
     * @{
     */

   //! Is there an aggretator of a particular type.
   bool isAggregator(const std::string & type) const;

   //! Is there an aggregator that uses a particular evaluation type
   template <typename EvalT>
   bool isAggregator() const;

   //! Is there an aggregator that is a particular type and uses a specified evaluation type
   template <typename EvalT>
   bool isAggregator(const std::string & type) const;

   /** @} */

   /** Get aggregator manager associated with a type string
     *
     * \note This method throws if no aggregator manager is found.
     */
   const AggregatorManager & getAggregatorManager(const std::string & type) const;

   Teuchos::RCP<LinearObjFactory<TraitsT> > getLinearObjFactory() const
   { return linObjFactory_; }

   Teuchos::RCP<UniqueGlobalIndexerBase > getGlobalIndexer() const
   { return globalIndexer_; }


   /** Add a set of panzer defined aggregators to the default aggregator manager.
     * In paricular the aggregator "Functional" is added as a Residual evaluator.
     */
   static void defineDefaultAggregators(ResponseAggregator_Manager<TraitsT> & aggMngr);
private:

   std::map<std::string,Teuchos::RCP<AggregatorManager> > aggregators_;

   Teuchos::RCP<UniqueGlobalIndexerBase > globalIndexer_;
   Teuchos::RCP<LinearObjFactory<TraitsT> > linObjFactory_;

   ResponseAggregator_Manager(const ResponseAggregator_Manager &);
};

} // end namespace panzer

#include "Panzer_ResponseAggregator_Manager_impl.hpp"

#endif
