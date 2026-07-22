// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_ROL_OBSERVERBASE_HPP
#define PIRO_ROL_OBSERVERBASE_HPP

#include <string>

#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

namespace Piro {

template <typename Scalar>
class ROL_ObserverBase {
public:

  virtual void observeSolution(
    double stamp, const Thyra::VectorBase<Scalar>& nonOverlappedSolution,
    const Teuchos::Ptr<const Thyra::MultiVectorBase<Scalar>>& nonOverlappedSolution_dxdp,
    const Teuchos::Ptr<const Thyra::VectorBase<Scalar>>& nonOverlappedSolutionDot,
    const Teuchos::Ptr<const Thyra::VectorBase<Scalar>>& nonOverlappedSolutionDotDot);

  virtual void observeSolution(
      double stamp, const Thyra::MultiVectorBase<Scalar>& nonOverlappedSolution, 
      const Teuchos::Ptr<const Thyra::MultiVectorBase<Scalar>>& nonOverlappedSolution_dxdp);

  virtual void parameterChanged(
      const std::string& param);

  virtual void parametersChanged();

  virtual void observeResponse(int iter);

  virtual ~ROL_ObserverBase() {}
};

template <typename Scalar>
void
ROL_ObserverBase<Scalar>::observeSolution(
    double stamp, const Thyra::VectorBase<Scalar>& /*nonOverlappedSolution*/,
    const Teuchos::Ptr<const Thyra::MultiVectorBase<Scalar>>& /*nonOverlappedSolution_dxdp*/,
    const Teuchos::Ptr<const Thyra::VectorBase<Scalar>>& /*nonOverlappedSolutionDot*/,
    const Teuchos::Ptr<const Thyra::VectorBase<Scalar>>& /*nonOverlappedSolutionDotDot*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ROL_ObserverBase<Scalar>::observeSolution(
      double stamp, const Thyra::MultiVectorBase<Scalar>& /*nonOverlappedSolution*/, 
      const Teuchos::Ptr<const Thyra::MultiVectorBase<Scalar>>& /*nonOverlappedSolution_dxdp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ROL_ObserverBase<Scalar>::parameterChanged(
      const std::string& param)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ROL_ObserverBase<Scalar>::parametersChanged()
{
  // Nothing to do by default
}

template <typename Scalar>
void
ROL_ObserverBase<Scalar>::observeResponse(
      int iter)
{
  // Nothing to do by default
}

} // namespace Piro

#endif /* PIRO_ROL_OBSERVERBASE_HPP */
