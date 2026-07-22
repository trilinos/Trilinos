// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCEDDYNAMICSTATIONARYCONTROLSOBJECTIVE_HPP
#define ROL_REDUCEDDYNAMICSTATIONARYCONTROLSOBJECTIVE_HPP

#include "ROL_Ptr.hpp"
#include "ROL_ReducedDynamicObjective.hpp"
#include "ROL_PartitionedVector.hpp"

namespace ROL
{

template <typename Real>
class ReducedDynamicStationaryControlsObjectiveHook
{
public:
  ReducedDynamicStationaryControlsObjectiveHook() {}
  virtual ~ReducedDynamicStationaryControlsObjectiveHook() {}
  virtual void
  preValue(const Vector<Real> &x) const = 0;
  virtual void
  postValue(Real val) const = 0;
  virtual void
  preGradient(const Vector<Real> &x) const = 0;
  virtual void
  postGradient(const Vector<Real> &g) const = 0;
};

/** @ingroup func_group
    \class ROL::ReducedDynamicStationaryControlsObjective
    \brief Defines the reduced time-dependent objective function interface
           for simulation-based optimization when the controls are stationary
           (i.e., not time-dependent).
*/
template <typename Real>
class ReducedDynamicStationaryControlsObjective : public Objective<Real>
{
public:
  using size_type = typename std::vector<Real>::size_type;

  ReducedDynamicStationaryControlsObjective(
      const Ptr<ReducedDynamicObjective<Real>> &red_dyn_obj,
      const Ptr<Vector<Real>> &x,
      const size_type Nt,
      const Ptr<ReducedDynamicStationaryControlsObjectiveHook<Real>> &hook = nullPtr)
      : red_dyn_obj_(red_dyn_obj), Nt_(Nt), hook_(hook)
  {
    x_dyn_ = PartitionedVector<Real>::create(*x, Nt_);
    v_dyn_ = PartitionedVector<Real>::create(*x, Nt_);
    g_dyn_ = PartitionedVector<Real>::create(x->dual(), Nt_);
  }

  virtual ~ReducedDynamicStationaryControlsObjective() {}

  void
  update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override
  {
    for (size_type i = 0; i < Nt_; ++i)
      (*x_dyn_)[i].set(x);
    red_dyn_obj_->update(*x_dyn_, true, iter);
  }  

  void
  update(const Vector<Real> &x, bool flag = true, int iter = -1) override
  {
    for (size_type i = 0; i < Nt_; ++i)
      (*x_dyn_)[i].set(x);
    red_dyn_obj_->update(*x_dyn_, flag, iter);
  }

  Real
  value(const Vector<Real> &x, Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->preValue(x);
    for (size_type i = 0; i < Nt_; ++i)
      (*x_dyn_)[i].set(x);
    Real val = red_dyn_obj_->value(*x_dyn_, tol);
    if (hook_ != nullPtr)
      hook_->postValue(val);
    return val;
  }

  void
  gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->preGradient(x);
    for (size_type i = 0; i < Nt_; ++i)
      (*x_dyn_)[i].set(x);
    red_dyn_obj_->gradient(*g_dyn_, *x_dyn_, tol);
    g.zero();
    for (size_type i = 0; i < Nt_; ++i)
      g.axpy(1.0, (*g_dyn_)[i]);
    if (hook_ != nullPtr)
      hook_->postGradient(g);
  }

  void
  hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override
  {
    for (size_type i = 0; i < Nt_; ++i)
    {
      (*x_dyn_)[i].set(x);
      (*v_dyn_)[i].set(v);
    }
    red_dyn_obj_->hessVec(*g_dyn_, *v_dyn_, *x_dyn_, tol);
    hv.zero();
    for (size_type i = 0; i < Nt_; ++i)
      hv.axpy(1.0, (*g_dyn_)[i]);
  }

private:
  Ptr<ReducedDynamicObjective<Real>> red_dyn_obj_;
  size_type Nt_;
  Ptr<ReducedDynamicStationaryControlsObjectiveHook<Real>> hook_;
  Ptr<PartitionedVector<Real>> x_dyn_;
  Ptr<PartitionedVector<Real>> v_dyn_;
  Ptr<PartitionedVector<Real>> g_dyn_;
};
} // namespace ROL

#endif
