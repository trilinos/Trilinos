// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Response_Probe_impl_hpp__
#define __Panzer_Response_Probe_impl_hpp__

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
Response_Probe<EvalT>::
Response_Probe(const std::string & responseName, MPI_Comm comm,
               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & linObjFact)
  : ResponseMESupport_Default<EvalT>(responseName,comm), value(0.0),
    have_probe(false), linObjFactory_(linObjFact)
{
  if(linObjFactory_!=Teuchos::null) {
    // requires thyra object factory
    thyraObjFactory_ = Teuchos::rcp_dynamic_cast<const panzer::ThyraObjFactory<double> >(linObjFactory_,true);
    setSolnVectorSpace(thyraObjFactory_->getThyraDomainSpace());

    // build a ghosted container, with a solution vector
    ghostedContainer_ = linObjFactory_->buildGhostedLinearObjContainer();

    // set ghosted container (work space for assembly)
    linObjFactory_->initializeGhostedContainer(panzer::LinearObjContainer::X,*ghostedContainer_);
  }
}

template <typename EvalT>
void Response_Probe<EvalT>::
initializeResponse()
{
  value = 0.0;
  have_probe = false;

  if(ghostedContainer_!=Teuchos::null) ghostedContainer_->initialize();
}

template <typename EvalT>
void Response_Probe<EvalT>::
scatterResponse()
{
  double glbValue = Sacado::scalarValue(value);

  // find the minimum processor who has the probe value
  int locProc = have_probe ? this->getComm()->getRank() : this->getComm()->getSize();
  int glbProc = 0;
  Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, Thyra::Ordinal(1), &locProc, &glbProc);

  TEUCHOS_ASSERT(glbProc < this->getComm()->getSize());

  // now broadcast the value from proc glbProc
  Teuchos::broadcast(*this->getComm(), glbProc, Thyra::Ordinal(1), &glbValue);

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
void Response_Probe<panzer::Traits::Jacobian>::
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
void Response_Probe<panzer::Traits::Hessian>::
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
void Response_Probe<panzer::Traits::Tangent>::
scatterResponse()
{
  const int n = value.size();
  const int num_deriv = this->numDeriv();
  TEUCHOS_ASSERT(n == 0 || n == num_deriv);
  if (n == 0)
    value.resize(num_deriv);

  // find the minimum processor who has the probe value
  if (num_deriv > 0) {
    int locProc = have_probe ? this->getComm()->getRank() : this->getComm()->getSize();
    int glbProc = 0;
    Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_MIN, Thyra::Ordinal(1), &locProc, &glbProc);

    TEUCHOS_ASSERT(glbProc < this->getComm()->getSize());

    // now broadcast the derivatives from proc glbProc
    Teuchos::broadcast(*this->getComm(), glbProc, Thyra::Ordinal(num_deriv), &value.fastAccessDx(0));
  }

  // copy data in vectors
#ifdef PANZER_HAVE_EPETRA_STACK
  if(this->useEpetra()) {
    // use epetra
    Epetra_MultiVector& deriv = this->getEpetraMultiVector();
    for (int i=0; i<num_deriv; ++i)
      deriv[i][0] = value.dx(i);
  }
  else
#endif
  {
    // use thyra
    TEUCHOS_ASSERT(this->useThyra());
    Thyra::ArrayRCP< Thyra::ArrayRCP<double> > deriv = this->getThyraMultiVector();
    for (int i=0; i<num_deriv; ++i)
      deriv[i][0] = value.dx(i);
  }
}

// Do nothing unless derivatives are actually required
template <typename EvalT>
void Response_Probe<EvalT>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & /* soln_vs */) { }

// derivatives are required for
template < >
void Response_Probe<panzer::Traits::Jacobian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
// derivatives are required for
template < >
void Response_Probe<panzer::Traits::Hessian>::
setSolnVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & soln_vs)
{
  setDerivativeVectorSpace(soln_vs);
}
#endif

// Do nothing unless derivatives are required
template <typename EvalT>
void Response_Probe<EvalT>::
adjustForDirichletConditions(const GlobalEvaluationData & /* localBCRows */, const GlobalEvaluationData & /* globalBCRows */) { }

// Do nothing unless derivatives are required
template < >
void Response_Probe<panzer::Traits::Jacobian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
// Do nothing unless derivatives are required
template < >
void Response_Probe<panzer::Traits::Hessian>::
adjustForDirichletConditions(const GlobalEvaluationData & localBCRows,const GlobalEvaluationData & globalBCRows)
{
  linObjFactory_->adjustForDirichletConditions(Teuchos::dyn_cast<const LinearObjContainer>(localBCRows),
                                               Teuchos::dyn_cast<const LinearObjContainer>(globalBCRows),
                                               *ghostedContainer_,true,true);
}
#endif

}

#endif
