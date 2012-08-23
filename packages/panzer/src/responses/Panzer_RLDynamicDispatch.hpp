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

#ifndef __Panzer_RLDynamicDispatch_hpp__
#define __Panzer_RLDynamicDispatch_hpp__

#include <string>

#include "Teuchos_Ptr.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Phalanx_TemplateManager.hpp"

namespace panzer {

// forward declaration
template <typename Traits> class ResponseLibrary;

/** \file Panzer_RLDynamicDispatch.hpp 
  * \brief Of use only by the internals of the response library.
  *
  * This provides a mechnism where a templated member function
  * can be called a run time using a dynamic call. This is specialized
  * to the ResponseLibrary, is meant to handle those functions that
  * use templating on the evaluation type.
  */

/** Abstract base class which acts as a surrogate for the
  * dynamic/static call.
  */
template <typename Traits>
class RLDynamicDispatchBase {
public:
   virtual ~RLDynamicDispatchBase() {}

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType) = 0;

   virtual const ResponseAggregatorBase<Traits> & getAggregator(const std::string & type,const std::string & evalType) const = 0;
};

/** A concrete instantiation of the response library
  * dynamic dispatch mechanism. This is the step that
  * calls the compile time polymorphic type from the dynamic string.
  * (CT stands for Compile Time)
  */
template <typename EvalT,typename Traits>
class RLDynamicDispatchCT : public RLDynamicDispatchBase<Traits> {
public:
   RLDynamicDispatchCT(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
      : responseLibrary_(rl) {}

   virtual ~RLDynamicDispatchCT() {}

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType)
   { 
      // sanity check
      TEUCHOS_TEST_FOR_EXCEPTION(PHX::TypeString<EvalT>::value!=evalType,std::logic_error,
                         "panzer::RLDynamicDispatch: Something is seriously wrong, dynamic evaluation type \""+evalType+"\" "
                         "is not equal to static type \""+PHX::TypeString<EvalT>::value+"\"!"); 
      responseLibrary_->template reserveVolumeResponse<EvalT>(rid,eBlock); 
   }

   virtual const ResponseAggregatorBase<Traits> & getAggregator(const std::string & type,const std::string & evalType) const
   { 
      // sanity check
      TEUCHOS_TEST_FOR_EXCEPTION(PHX::TypeString<EvalT>::value!=evalType,std::logic_error,
                         "panzer::RLDynamicDispatch: Something is seriously wrong, dynamic evaluation type \""+evalType+"\" "
                         "is not equal to static type \""+PHX::TypeString<EvalT>::value+"\"!"); 
      return responseLibrary_->template getAggregator<EvalT>(type);
   }

private:
   RLDynamicDispatchCT();
   RLDynamicDispatchCT(const RLDynamicDispatchCT<EvalT,Traits> &);

   Teuchos::Ptr<ResponseLibrary<Traits> > responseLibrary_;
};

template <typename Traits>
class RLDynamicDispatch : public RLDynamicDispatchBase<Traits> {
public:
   RLDynamicDispatch();

   void buildObjects(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl);

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType)
   { dynToStaticMap_[evalType]->reserveVolumeResponse(rid,eBlock,evalType); }

   virtual const ResponseAggregatorBase<Traits> & getAggregator(const std::string & type,const std::string & evalType) const
   { return (dynToStaticMap_.find(evalType)->second)->getAggregator(type,evalType); }

   const std::vector<std::string> & getEvaluationNames() const
   { return evaluationNames_; }

private:
   /** Class to turn an mpl::vector into a string vector
     * using the <code>PHX::TypeString</code> class.
     */
   struct NameFill {
      //! mutable so that we can modify this and use the for_each in Sacado
      std::vector<std::string> & vec;

      NameFill(std::vector<std::string> & v) : vec(v) {}
      NameFill(const NameFill & n) : vec(n.vec) {}

      template <typename T> void operator()(T x) const
      { vec.push_back(PHX::TypeString<T>::value); }
      
   private:
      NameFill();
   };

   //! For use with a template manager for building dynamic dispatch CT objects
   struct RLDynDispCTBuilder {
      Teuchos::Ptr<ResponseLibrary<Traits> > responseLibrary_;

      RLDynDispCTBuilder(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
         : responseLibrary_(rl) {}
 
      template <typename EvalT>
      Teuchos::RCP<RLDynamicDispatchBase<Traits> > build() const
      { return Teuchos::rcp(new RLDynamicDispatchCT<EvalT,Traits>(responseLibrary_)); }
   };

   std::vector<std::string> evaluationNames_;

   std::map<std::string,Teuchos::RCP<RLDynamicDispatchBase<Traits> > > dynToStaticMap_;
};

template <typename Traits>
RLDynamicDispatch<Traits>::RLDynamicDispatch()
{
   typedef typename Traits::EvalTypes EvalTypes;
   Sacado::mpl::for_each<EvalTypes>(NameFill(evaluationNames_));
}

template <typename Traits>
void RLDynamicDispatch<Traits>::buildObjects(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
{
   typedef PHX::TemplateManager<typename Traits::EvalTypes,RLDynamicDispatchBase<Traits>,RLDynamicDispatchCT<_,Traits> > 
           DynDispatchTM;

   DynDispatchTM dynDispatchTM;
   dynDispatchTM.buildObjects(RLDynDispCTBuilder(rl));

   // build map from string to dynamic dispatch object
   std::vector<std::string>::const_iterator nameItr = evaluationNames_.begin();
   typename DynDispatchTM::iterator dynDspItr(dynDispatchTM.begin());
   for(; nameItr!=evaluationNames_.end() && dynDspItr!=dynDispatchTM.end();
       ++nameItr,++dynDspItr) {
      dynToStaticMap_[*nameItr] = dynDspItr.rcp();
   }
}

}

#endif
