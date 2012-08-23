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

#ifndef __Panzer_ResponseContainer_hpp__
#define __Panzer_ResponseContainer_hpp__

#include <vector>
#include <set>
#include <string>
#include <list>

#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"

#include "Panzer_Response.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_PhysicsBlock.hpp"

namespace panzer {

/** Base class for use with the Phalanx TemplateManager.
  * This object will serve as an accessor to the various
  * response objects.
  */
template <typename TraitsT>
class ResponseContainerBase {
public:
   ResponseContainerBase() {}

   ResponseContainerBase(const Teuchos::RCP<const ResponseLibrary<TraitsT> > & rl)
      : responseLibrary_(rl) {}

   void setResponseLibrary(const Teuchos::RCP<const ResponseLibrary<TraitsT> > & rl)
   { responseLibrary_ = rl; }

   virtual ~ResponseContainerBase() {}

   /** \defgroup Methods common to all ResponseContainers
     * These methods are implementations of a basic set of functionality
     * common to all ResponseContainer classes.
     * @{
     */

   //! Reserve a particular response
   void reserve(const ResponseId & rid)
   { responses_.insert(rid); }

   //! Does this container have a response of a specified name
   bool contains(const ResponseId & rid) const
   {
      std::set<ResponseId>::const_iterator itr = std::find(responses_.begin(),responses_.end(),rid);
      return itr!=responses_.end();
   }

   //! get a list of all reserved responses
   const std::set<ResponseId> & getReserved() const
   { return responses_; }

   //! Clear response containers data
   void clear()
   {
      typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator itr;
      for(itr=responseDataObjs_.begin();itr!=responseDataObjs_.end();++itr)
         itr->second->reinitializeData();
   }

   // provide reference only access to the response library
   const ResponseLibrary<TraitsT> & getResponseLibrary() const 
   { return *responseLibrary_; }

   Teuchos::RCP<ResponseData<TraitsT> > getResponseData(const std::string & type) const
   {
      typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::const_iterator itr
           = responseDataObjs_.find(type);
      if(itr==endDataObjs())
         return Teuchos::null;
      return itr->second;
   }

   void print(std::ostream & os) const 
   {
      const std::string & evalType = getResponseEvalType();
      const std::set<ResponseId> & reserved = getReserved();

      os << "Eval Type = " << evalType << ", Responses = ";
      for(std::set<ResponseId>::const_iterator itr=reserved.begin();
          itr!=reserved.end();++itr) {
         os << *itr << ", ";
      }
   }

   /** @} */

   /** \defgroup Iterator access to the response data.
     * @{
     */

   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator
   beginDataObjs() { return responseDataObjs_.begin(); }

   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator
   endDataObjs()   { return responseDataObjs_.end(); }

   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::const_iterator
   beginDataObjs() const { return responseDataObjs_.begin(); } 

   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::const_iterator
   endDataObjs() const   { return responseDataObjs_.end(); }

   /** @} */

   /** \defgroup Pure virtual funcions for setting up field managers
     * @{
     */

   //! What is the evaluation type for this response
   virtual const std::string & getResponseEvalType() const = 0;

   /** Register responses with a particular field manager.
     * This essentially builds a ScatterResponse object 
     * and registers it as a required field in the field manager.
     */
   virtual void registerResponses(PHX::FieldManager<TraitsT> & fm,const PhysicsBlock & pb,const Teuchos::ParameterList & pl) = 0;

   //! Using reserved data build the response data objects
   virtual void buildResponseDataObjs() = 0;

   //! perform global reduction on this set of response data
   virtual void globalReduction(const Teuchos::Comm<int> & comm) = 0;

   /** @} */

protected:

   /** Insert a response data object into the container
     *
     * \note This method wipes out any previously data object.
     */
   void insertResponseDataObj(const std::string & type,const Teuchos::RCP<ResponseData<TraitsT> > & rd)
   { responseDataObjs_[type] = rd; }

private:
   //! Response library for accesing the aggregators
   Teuchos::RCP<const ResponseLibrary<TraitsT> > responseLibrary_;

   //! Responses contained in this container.
   std::set<ResponseId> responses_;

   //! Stores all response data objects.
   std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > > responseDataObjs_; 

   // hide undesirable constructors
   ResponseContainerBase(const ResponseContainerBase<TraitsT> &);
};

template <typename TraitsT>
std::ostream & operator<<(std::ostream & os,const ResponseContainerBase<TraitsT> & container)
{ 
   container.print(os); 
   return os; 
}

// forward declaration
template <typename EvalT,typename TraitsT> 
class ResponseContainer 
   : public ResponseContainerBase<TraitsT> {

public:
   ResponseContainer()
      : ResponseContainerBase<TraitsT>(), responseDataObjsBuilt_(false) {}

   ResponseContainer(const Teuchos::RCP<const ResponseLibrary<TraitsT> > & rl)
      : ResponseContainerBase<TraitsT>(rl), responseDataObjsBuilt_(false) {}

   //! What is the evaluation type for this response
   virtual const std::string & getResponseEvalType() const 
   { return PHX::TypeString<EvalT>::value; }

   /** Register responses with a particular field manager.
     * This essentially builds a ScatterResponse object 
     * and registers it as a required field in the field manager.
     */
   virtual void registerResponses(PHX::FieldManager<TraitsT> & fm,const PhysicsBlock & pb,const Teuchos::ParameterList & pl);

   //! Using reserved data build the response data objects
   void buildResponseDataObjs();

   //! perform global reduction on this set of response data
   virtual void globalReduction(const Teuchos::Comm<int> & comm);

private:
   bool responseDataObjsBuilt_;
};

template <typename EvalT,typename TraitsT> 
void ResponseContainer<EvalT,TraitsT>::
registerResponses(PHX::FieldManager<TraitsT> & fm,const PhysicsBlock & pb,const Teuchos::ParameterList & pl)
{
   using Teuchos::RCP;

   if(!responseDataObjsBuilt_) 
      buildResponseDataObjs();

   const ResponseLibrary<TraitsT> & rLibrary = this->getResponseLibrary();

   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator beg
      = this->beginDataObjs();
   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator end
      = this->endDataObjs();

   for(typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator itr=beg;
       itr!=end;++itr) {
      // extract pieces
      std::string type = itr->first;
      RCP<ResponseData<TraitsT> > dataObj = itr->second;

      // use aggregator register required evaluators
      const ResponseAggregatorBase<Traits> & aggregator 
            = rLibrary.template getAggregator<EvalT>(type);
      aggregator.registerAndRequireEvaluators(fm,dataObj,pb,pl);
   }
}

template <typename EvalT,typename TraitsT> 
void ResponseContainer<EvalT,TraitsT>::
buildResponseDataObjs()
{
   const std::set<ResponseId> & reserved = this->getReserved();
   std::map<std::string,std::vector<std::string> > fields;

   const ResponseLibrary<TraitsT> & rLibrary = this->getResponseLibrary();
 
   for(std::set<ResponseId>::const_iterator itr=reserved.begin();
       itr!=reserved.end();++itr) {
      const ResponseId & response = *itr;

      // make sure nothing wacky is going on
      TEUCHOS_ASSERT(rLibrary.template isResponseType<EvalT>(response.type));

      // add in the response
      fields[response.type].push_back(response.name);
   }

   // build data objects for each type
   for(std::map<std::string,std::vector<std::string> >::const_iterator itr=fields.begin();
       itr!=fields.end();++itr) {

      std::string type = itr->first;

      // this should never fail because of the assertion above
      const ResponseAggregatorBase<Traits> & aggregator 
            = rLibrary.template getAggregator<EvalT>(type);

      // add a new response data object to this container associated with the type
      this->insertResponseDataObj(type,aggregator.buildResponseData(itr->second));
   }

   responseDataObjsBuilt_ = true;
}

//! perform global reduction on this set of response data
template <typename EvalT,typename TraitsT> 
void ResponseContainer<EvalT,TraitsT>::
globalReduction(const Teuchos::Comm<int> & comm)
{
   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator beg
      = this->beginDataObjs();
   typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator end
      = this->endDataObjs();

   const ResponseLibrary<TraitsT> & rLibrary = this->getResponseLibrary();

   for(typename std::map<std::string,Teuchos::RCP<ResponseData<TraitsT> > >::iterator itr=beg;
       itr!=end;++itr) {
      // extract pieces
      std::string type = itr->first;
      Teuchos::RCP<ResponseData<TraitsT> > dataObj = itr->second;

      // use aggregator register required evaluators
      const ResponseAggregatorBase<Traits> & aggregator 
            = rLibrary.template getAggregator<EvalT>(type);
      aggregator.globalReduction(comm,*dataObj);
   }
}

}

#endif
