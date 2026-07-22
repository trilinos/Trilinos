// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Response_Functional_impl_hpp__
#define __Panzer_Response_Functional_impl_hpp__

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Epetra_LocalMap.h"
#endif

#include "Sacado_Traits.hpp"

namespace panzer {

template <typename EvalT>
void Response_Functional<EvalT>::
scatterResponse()
{
  double locValue = Sacado::scalarValue(value);
  double glbValue = 0.0;

  // do global summation
  Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_SUM, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);

  value = glbValue;

  // built data in vectors
#ifdef PANZER_HAVE_EPETRA_STACK
  if(this->useEpetra()) {
    // use epetra
    this->getEpetraVector()[0] = glbValue;
  }
  else
 #endif
  {
    // use thyra
    TEUCHOS_ASSERT(this->useThyra());

    this->getThyraVector()[0] = glbValue;
  }
}

template < >
void Response_Functional<panzer::Traits::Jacobian>::
scatterResponse()
{
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::MultiVectorBase<double> > dgdx_unique = getDerivative();

  // if its null, don't do anything
  if(dgdx_unique==Teuchos::null)
    return;

  uniqueContainer_ = linObjFactory_->buildLinearObjContainer();
  Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(uniqueContainer_)->set_x_th(dgdx_unique->col(0));

  linObjFactory_->ghostToGlobalContainer(*ghostedContainer_,*uniqueContainer_,LinearObjContainer::X);

  uniqueContainer_ = Teuchos::null;
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
void Response_Functional<panzer::Traits::Hessian>::
scatterResponse()
{
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::MultiVectorBase<double> > dgdx_unique = getDerivative();

  // if its null, don't do anything
  if(dgdx_unique==Teuchos::null)
    return;

  uniqueContainer_ = linObjFactory_->buildLinearObjContainer();
  Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(uniqueContainer_)->set_x_th(dgdx_unique->col(0));

  linObjFactory_->ghostToGlobalContainer(*ghostedContainer_,*uniqueContainer_,LinearObjContainer::X);

  uniqueContainer_ = Teuchos::null;
}
#endif

template < >
void Response_Functional<panzer::Traits::Tangent>::
scatterResponse()
{
  const int n = value.size();
  const int num_deriv = this->numDeriv();
  TEUCHOS_ASSERT(n == 0 || n == num_deriv);
  ScalarT glbValue = ScalarT(num_deriv, 0.0);

  // do global summation -- it is possible to do the reduceAll() on the Fad's directly, but it is somewhat
  // complicated for DFad (due to temporaries that might get created).  Since this is just a sum, it is
  // easier to do the reduction for each value and derivative component.
  Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_SUM, Thyra::Ordinal(1), &value.val(), &glbValue.val());
  if (num_deriv > 0)
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_SUM, Thyra::Ordinal(n), value.dx(),  &glbValue.fastAccessDx(0));

  value = glbValue;

  // copy data in vectors
#ifdef PANZER_HAVE_EPETRA_STACK
  if(this->useEpetra()) {
    // use epetra
    Epetra_MultiVector& deriv = this->getEpetraMultiVector();
    for (int i=0; i<num_deriv; ++i)
      deriv[i][0] = glbValue.dx(i);
  }
  else
#endif
  {
    // use thyra
    TEUCHOS_ASSERT(this->useThyra());
    Thyra::ArrayRCP< Thyra::ArrayRCP<double> > deriv = this->getThyraMultiVector();
    for (int i=0; i<num_deriv; ++i)
      deriv[i][0] = glbValue.dx(i);
  }
}

// Do nothing unless derivatives are actually required
template <typename EvalT>
void Response_Functional<EvalT>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & /* soln_vs */) { }

// derivatives are required for
template < >
void Response_Functional<panzer::Traits::Jacobian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}

// derivatives are required for
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
void Response_Functional<panzer::Traits::Hessian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}
#endif

// Do nothing unless derivatives are required
template <typename EvalT>
void Response_Functional<EvalT>::
adjustForDirichletConditions(const GlobalEvaluationData & /* localBCRows */, const GlobalEvaluationData & /* globalBCRows */) { }

// Do nothing unless derivatives are required
template < >
void Response_Functional<panzer::Traits::Jacobian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
// Do nothing unless derivatives are required
template < >
void Response_Functional<panzer::Traits::Hessian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}
#endif

}

#endif
