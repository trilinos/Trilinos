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

#ifndef __Panzer_Response_Functional_impl_hpp__
#define __Panzer_Response_Functional_impl_hpp__

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Epetra_LocalMap.h"

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
  if(this->useEpetra()) {
    // use epetra
    Epetra_MultiVector& deriv = this->getEpetraMultiVector();
    for (int i=0; i<num_deriv; ++i)
      deriv[i][0] = glbValue.dx(i);
  }
  else {
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
