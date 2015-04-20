#ifndef __Panzer_ResponseMESupportBuilderBase_hpp__
#define __Panzer_ResponseMESupportBuilderBase_hpp__

#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"

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
  virtual void setDerivativeInformationBase(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory,
                                            const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & globalIndexer) = 0;

  /** Using a panzer::Residual evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildValueFactory() const = 0;

  /** Using a panzer::Jacobian evaluation type build the REFB for this
    * response.
    */
  virtual Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> buildDerivativeFactory() const = 0;

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

}

#endif
