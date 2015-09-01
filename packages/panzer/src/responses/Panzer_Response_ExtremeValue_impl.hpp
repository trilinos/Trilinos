#ifndef __Panzer_Response_ExtremeValue_impl_hpp__
#define __Panzer_Response_ExtremeValue_impl_hpp__

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Epetra_LocalMap.h"

#include "Sacado_Traits.hpp"

namespace panzer {

template <typename EvalT>
void Response_ExtremeValue<EvalT>::
scatterResponse() 
{
  double locValue = Sacado::ScalarValue<ScalarT>::eval(value);
  double glbValue = 0.0;

  // do global summation
  if(useMax_)
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MAX, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);
  else
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);

  value = glbValue;

  // built data in vectors
  if(this->useEpetra()) {
    // use epetra 
    this->getEpetraVector()[0] = glbValue;
  }
  else {
    // use thyra
    TEUCHOS_ASSERT(this->useThyra());

    this->getThyraVector()[0] = glbValue;
  }
}

template < >
void Response_ExtremeValue<panzer::Traits::Jacobian>::
scatterResponse() 
{
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::MultiVectorBase<double> > dgdx_unique = getDerivative();

  uniqueContainer_ = linObjFactory_->buildLinearObjContainer();
  Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(uniqueContainer_)->set_x_th(dgdx_unique->col(0));

  linObjFactory_->ghostToGlobalContainer(*ghostedContainer_,*uniqueContainer_,LinearObjContainer::X);

  uniqueContainer_ = Teuchos::null;
}

// Do nothing unless derivatives are actually required
template <typename EvalT>
void Response_ExtremeValue<EvalT>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs) { }

// derivatives are required for 
template < >
void Response_ExtremeValue<panzer::Traits::Jacobian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs) 
{ 
  setDerivativeVectorSpace(soln_vs);
}

// Do nothing unless derivatives are required
template <typename EvalT>
void Response_ExtremeValue<EvalT>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows) { }

// Do nothing unless derivatives are required
template < >
void Response_ExtremeValue<panzer::Traits::Jacobian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows) 
{ 
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}

}

#endif
