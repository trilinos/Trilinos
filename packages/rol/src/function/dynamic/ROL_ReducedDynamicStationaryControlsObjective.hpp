// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
