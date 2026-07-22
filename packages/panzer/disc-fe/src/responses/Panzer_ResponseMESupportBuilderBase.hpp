// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
