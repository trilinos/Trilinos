// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_SOLUTIONOBSERVERBASE_HPP
#define PIRO_SOLUTIONOBSERVERBASE_HPP

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

#include "Thyra_VectorBase.hpp"

namespace Piro {

template <typename Scalar, typename VectorType>
class SolutionObserverBase {
public:
  virtual void observeResponse(
      int j,
      const Teuchos::RCP<Thyra::ModelEvaluatorBase::OutArgs<Scalar> >& outArgs,
      const Teuchos::RCP<Teuchos::Array<Teuchos::RCP<VectorType> > > &responses,
      const Teuchos::RCP<VectorType> &g);

  virtual ~SolutionObserverBase() {}
};

template <typename Scalar, typename VectorType>
void
SolutionObserverBase<Scalar, VectorType>::observeResponse(
      int /* j */,
      const Teuchos::RCP<Thyra::ModelEvaluatorBase::OutArgs<Scalar> >& /* outArgs */,
      const Teuchos::RCP<Teuchos::Array<Teuchos::RCP<VectorType> > > &/* responses */,
      const Teuchos::RCP<VectorType> &/* g */)
{
  // Nothing to do by default
}

} // namespace Piro

#endif /* PIRO_SOLUTIONOBSERVERBASE_HPP */
