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

#ifndef __Panzer_ResponseMESupportBuilderBase_hpp__
#define __Panzer_ResponseMESupportBuilderBase_hpp__

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"

namespace panzer {

/** This class is used by the model evaluator and it supports setting up derivative information for a 
  * response. In particular, it provides a mechanism for defining which distributed parameters are used for
  * compute derivatives.
  */
class ResponseMESupportBuilderBase {
public:
  virtual ~ResponseMESupportBuilderBase() {}
 
  /** This method controls how the derivative vector is allocated and 
    * scattered. The idea here is a Response can have different partial
    * derivatives and this provides the mechanism for supporting that.
    */
  virtual void setDerivativeInformation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory) = 0;

  /** Using a panzer::Residual evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildValueFactory() const = 0;

  /** Using a panzer::Jacobian evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildDerivativeFactory() const = 0;

  /** Using a panzer::Tangent evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildTangentFactory() const {
    return Teuchos::null;
  }

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  /** Using a panzer::Tangent evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildHessianFactory() const = 0;
#endif

  /** Satisfy the required interface for the builder used in the "addResponse" function
    * in the ResponseLibrary.
    */
  template <typename T>
  inline Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::null; }
};

template < >
inline Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> ResponseMESupportBuilderBase::build<panzer::Traits::Residual>() const
{ return buildValueFactory(); }

template < >
inline Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> ResponseMESupportBuilderBase::build<panzer::Traits::Jacobian>() const
{ return buildDerivativeFactory(); }

template < >
inline Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> ResponseMESupportBuilderBase::build<panzer::Traits::Tangent>() const
{ return buildTangentFactory(); }

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
inline Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> ResponseMESupportBuilderBase::build<panzer::Traits::Hessian>() const
{ return buildHessianFactory(); }
#endif

}

#endif
