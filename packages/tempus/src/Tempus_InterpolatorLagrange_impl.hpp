//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_InterpolatorLagrange_impl_hpp
#define Tempus_InterpolatorLagrange_impl_hpp

#include <algorithm>
#include <iterator>

#include "Teuchos_ScalarTraits.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <class Scalar>
void InterpolatorLagrange<Scalar>::interpolate(
    const Scalar& t, SolutionState<Scalar>* state_out) const
{
  int n = (*nodes_).size();
  TEUCHOS_ASSERT(n > 0);
  if (n >= order_ + 1)
    lagrange(order_, t, state_out);
  else
    lagrange(n - 1, t, state_out);
}

template <class Scalar>
void InterpolatorLagrange<Scalar>::setParameterList(
    const Teuchos::RCP<Teuchos::ParameterList>& paramList)
{
  pl_ = Teuchos::parameterList();
  if (paramList == Teuchos::null)
    *pl_ = *(this->getValidParameters());
  else
    *pl_ = *paramList;
  pl_->validateParametersAndSetDefaults(*this->getValidParameters());
  order_ = pl_->get<int>("Order");
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
InterpolatorLagrange<Scalar>::getNonconstParameterList()
{
  return pl_;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
InterpolatorLagrange<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> tmp = pl_;
  pl_                                      = Teuchos::null;
  return tmp;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
InterpolatorLagrange<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> tmp = Teuchos::parameterList();
  tmp->set<std::string>("Interpolator Type", "Lagrange");
  tmp->set<int>("Order", 0);
  return tmp;
}

template <class Scalar>
void InterpolatorLagrange<Scalar>::lagrange(
    const int p, const Scalar& t, SolutionState<Scalar>* state_out) const
{
  using Teuchos::is_null;
  using Teuchos::RCP;
  using Thyra::V_StVpStV;
  using Thyra::VectorBase;

  // Here we assume we have at least p nodes
  int n = nodes_->size();
  TEUCHOS_ASSERT(n >= p);

  const Scalar t_begin = (*nodes_)[0]->getTime();
  const Scalar t_final = (*nodes_)[n - 1]->getTime();

  // Find largest index i such that
  //      (*nodes)[i]->getTime() <= t
  int i;
  if (t <= t_begin)
    i = 0;
  else if (t >= t_final)
    i = n - 1;
  else {
    auto it = std::find_if(nodes_->begin(), nodes_->end(),
                           [=](const RCP<SolutionState<Scalar> >& s) {
                             return t <= s->getTime();
                           });
    i       = std::distance(nodes_->begin(), it) - 1;
  }
  TEUCHOS_ASSERT(i >= 0 && i <= n - 1);

  // First we check for exact node matches:
  if (floating_compare_equals((*nodes_)[i]->getTime(), t, t_final)) {
    state_out->copy((*nodes_)[i]);
    return;
  }
  if (i < n - 1 &&
      floating_compare_equals((*nodes_)[i + 1]->getTime(), t, t_final)) {
    state_out->copy((*nodes_)[i + 1]);
    return;
  }

  // Put t roughly in the middle of the interpolation window
  int k = i - p / 2;
  if (k < 0) k = 0;
  if (k + p + 1 > n) k = n - p - 1;
  TEUCHOS_ASSERT(k >= 0 && k + p + 1 <= n);

  RCP<VectorBase<Scalar> > x       = state_out->getX();
  RCP<VectorBase<Scalar> > xdot    = state_out->getXDot();
  RCP<VectorBase<Scalar> > xdotdot = state_out->getXDotDot();

  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  Thyra::assign(x.ptr(), zero);
  if (!is_null(xdot)) Thyra::assign(xdot.ptr(), zero);
  if (!is_null(xdotdot)) Thyra::assign(xdotdot.ptr(), zero);

  for (int j = 0; j <= p; ++j) {
    const Scalar tj                         = (*nodes_)[k + j]->getTime();
    RCP<const VectorBase<Scalar> > xj       = (*nodes_)[k + j]->getX();
    RCP<const VectorBase<Scalar> > xdotj    = (*nodes_)[k + j]->getXDot();
    RCP<const VectorBase<Scalar> > xdotdotj = (*nodes_)[k + j]->getXDotDot();

    Scalar num = 1.0;
    Scalar den = 1.0;
    for (int l = 0; l <= p; ++l) {
      if (l != j) {
        const Scalar tl = (*nodes_)[k + l]->getTime();
        num *= t - tl;
        den *= tj - tl;
      }
    }
    const Scalar a = num / den;
    Thyra::Vp_StV(x.ptr(), a, *xj);
    if (!is_null(xdot)) Thyra::Vp_StV(xdot.ptr(), a, *xdotj);
    if (!is_null(xdotdot)) Thyra::Vp_StV(xdotdot.ptr(), a, *xdotdotj);
  }

  // Set meta-data for interpolated state
  state_out->getMetaData()->copy((*nodes_)[i]->getMetaData());
  state_out->getMetaData()->setTime(t);
  state_out->getMetaData()->setDt(t - (*nodes_)[i]->getTime());
  state_out->getMetaData()->setIsInterpolated(true);
  // What else should we set??
}

}  // namespace Tempus

#endif  // Tempus_InterpolatorLagrange_impl_hpp
