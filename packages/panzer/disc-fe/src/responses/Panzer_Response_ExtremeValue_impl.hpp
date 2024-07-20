// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Response_ExtremeValue_impl_hpp__
#define __Panzer_Response_ExtremeValue_impl_hpp__

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
void Response_ExtremeValue<EvalT>::
scatterResponse()
{
  double locValue = Sacado::scalarValue(value);
  double glbValue = 0.0;

  // do global summation
  if(useMax_)
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MAX, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);
  else
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);

  value = glbValue;

#ifdef PANZER_HAVE_EPETRA_STACK
  // built data in vectors
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

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
void Response_ExtremeValue<panzer::Traits::Hessian>::
scatterResponse()
{
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::MultiVectorBase<double> > dgdx_unique = getDerivative();

  uniqueContainer_ = linObjFactory_->buildLinearObjContainer();
  Teuchos::rcp_dynamic_cast<ThyraObjContainer<double> >(uniqueContainer_)->set_x_th(dgdx_unique->col(0));

  linObjFactory_->ghostToGlobalContainer(*ghostedContainer_,*uniqueContainer_,LinearObjContainer::X);

  uniqueContainer_ = Teuchos::null;
}
#endif

template < >
void Response_ExtremeValue<panzer::Traits::Tangent>::
scatterResponse()
{
  const int n = value.size();
  const int num_deriv = this->numDeriv();
  TEUCHOS_ASSERT(n == 0 || n == num_deriv);
  ScalarT glbValue = ScalarT(num_deriv, 0.0);

  // do global min/max on value
  if(useMax_)
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MAX, Thyra::Ordinal(1), &value.val(), &glbValue.val());
  else
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, Thyra::Ordinal(1), &value.val(), &glbValue.val());

  // find the minimum processor who's local value == the global min/max value
  if (num_deriv > 0) {
    int locProc = value.val() == glbValue.val() ? this->getComm()->getRank() : this->getComm()->getSize();
    int glbProc = 0;
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, Thyra::Ordinal(1), &locProc, &glbProc);

    // now broadcast the derivatives from proc glbProc
    Teuchos::broadcast(*this->getComm(), glbProc, Thyra::Ordinal(n), &glbValue.fastAccessDx(0));
  }

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
void Response_ExtremeValue<EvalT>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & /* soln_vs */) { }

// derivatives are required for
template < >
void Response_ExtremeValue<panzer::Traits::Jacobian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
// derivatives are required for
template < >
void Response_ExtremeValue<panzer::Traits::Hessian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}
#endif

// Do nothing unless derivatives are required
template <typename EvalT>
void Response_ExtremeValue<EvalT>::
adjustForDirichletConditions(const GlobalEvaluationData & /* localBCRows */, const GlobalEvaluationData & /* globalBCRows */) { }

// Do nothing unless derivatives are required
template < >
void Response_ExtremeValue<panzer::Traits::Jacobian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
// Do nothing unless derivatives are required
template < >
void Response_ExtremeValue<panzer::Traits::Hessian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}
#endif

}

#endif
