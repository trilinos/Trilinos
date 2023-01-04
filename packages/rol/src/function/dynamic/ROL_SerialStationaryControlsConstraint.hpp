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

#ifndef ROL_SERIALSTATIONARYCONTROLSCONSTRAINT_HPP
#define ROL_SERIALSTATIONARYCONTROLSCONSTRAINT_HPP

#include "ROL_Ptr.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_PartitionedVector.hpp"

namespace ROL
{

template <typename Real>
class SerialStationaryControlsConstraintHook
{
public:
  SerialStationaryControlsConstraintHook() {}
  virtual ~SerialStationaryControlsConstraintHook() {}
  virtual void
  trace_function(const std::string &name, const Vector<Real> &z) const = 0;
};

/** @ingroup func_group
    \class ROL::SerialStationaryControlsConstraint
    \brief Wrapper for SerialConstraint for when the controls are stationary
           (i.e., not time-dependent).
*/
template <typename Real>
class SerialStationaryControlsConstraint : public Constraint_SimOpt<Real>
{
public:
  using size_type = typename std::vector<Real>::size_type;

  SerialStationaryControlsConstraint(
      const Ptr<SerialConstraint<Real>> &serial_con,
      const Ptr<Vector<Real>> &z,
      const size_type Nt,
      const Ptr<SerialStationaryControlsConstraintHook<Real>> &hook = nullPtr)
      : serial_con_(serial_con), Nt_(Nt), hook_(hook)
  {
    z_dyn_ = PartitionedVector<Real>::create(*z, Nt_);
    v_dyn_ = PartitionedVector<Real>::create(*z, Nt_);
    g_dyn_ = PartitionedVector<Real>::create(z->dual(), Nt_);
  }

  virtual ~SerialStationaryControlsConstraint() {}

  void
  solve(Vector<Real> &c, Vector<Real> &u, const Vector<Real> &z, Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("solve", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->solve(c, u, *z_dyn_, tol);
  }

  void
  update(const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("value", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->update(u, *z_dyn_, flag, iter);
  }

  void
  value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("value", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->value(c, u, *z_dyn_, tol);
  }

  void
  applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
      override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyJacobian_1", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyJacobian_1(jv, v, u, *z_dyn_, tol);
  }

  void
  applyInverseJacobian_1(
      Vector<Real> &ijv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyInverseJacobian_1", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyInverseJacobian_1(ijv, v, u, *z_dyn_, tol);
  }

  // Serial constraint doesn't overload this one, so it is hidden
  // void applyAdjointJacobian_1(Vector<Real>& ajv,
  //                             const Vector<Real>& v,
  //                             const Vector<Real>& u,
  //                             const Vector<Real>& z,
  //                             Real& tol) override
  // {
  //   for (size_type i=0; i<Nt_; ++i)
  //     (*z_dyn_)[i].set(z);
  //   serial_con_->applyAdjointJacobian_1(ajv, v, u, *z_dyn_, tol);
  // }

  void
  applyAdjointJacobian_1(
      Vector<Real> &ajv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      const Vector<Real> &dualv,
      Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyAdjointJacobian_1", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyAdjointJacobian_1(ajv, v, u, *z_dyn_, dualv, tol);
  }

  void
  applyInverseAdjointJacobian_1(
      Vector<Real> &iajv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyInverseAdjointJacobian_1", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyInverseAdjointJacobian_1(iajv, v, u, *z_dyn_, tol);
  }

  void
  applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
      override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyJacobian_2", z);
    for (size_type i = 0; i < Nt_; ++i)
    {
      (*z_dyn_)[i].set(z);
      (*v_dyn_)[i].set(v);
    }
    serial_con_->applyJacobian_2(jv, *v_dyn_, u, *z_dyn_, tol);
  }

  void
  applyAdjointJacobian_2(
      Vector<Real> &ajv,
      const Vector<Real> &v,
      const Vector<Real> &u,
      const Vector<Real> &z,
      Real &tol) override
  {
    if (hook_ != nullPtr)
      hook_->trace_function("applyAdjointJacobian_2", z);
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyAdjointJacobian_2(*g_dyn_, v, u, *z_dyn_, tol);
    ajv.zero();
    for (size_type i = 0; i < Nt_; ++i)
      ajv.axpy(1.0, (*g_dyn_)[i]);
  }

  // Serial constraint doesn't overload this one, so it is hidden
  // void applyAdjointJacobian_2(Vector<Real>& ajv,
  //                             const Vector<Real>& v,
  //                             const Vector<Real>& u,
  //                             const Vector<Real>& z,
  //                             const Vector<Real>& dualv,
  //                             Real& tol) override
  // {
  //   for (size_type i=0; i<Nt_; ++i)
  //     (*z_dyn_)[i].set(z);
  //   serial_con_->applyAdjointJacobian_2(*g_dyn_, v, u, *z_dyn_, dualv, tol);
  //   ajv.zero();
  //   for (size_type i=0; i<Nt_; ++i)
  //     ajv.axpy(1.0, (*g_dyn_)[i]);
  // }

  /* The following I believe are correct, but are untested
  void applyAdjointHessian_11(Vector<Real>& ahwv,
                              const Vector<Real>& w,
                              const Vector<Real>& v,
                              const Vector<Real>& u,
                              const Vector<Real>& z,
                              Real& tol) override
  {
    for (size_type i=0; i<Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyAdjointHessian_11(ahwv, w, v, u, *z_dyn_, tol);
  }

  void applyAdjointHessian_12(Vector<Real>& ahwv,
                              const Vector<Real>& w,
                              const Vector<Real>& v,
                              const Vector<Real>& u,
                              const Vector<Real>& z,
                              Real& tol) override
  {
    for (size_type i=0; i<Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_con_->applyAdjointHessian_12(*g_dyn_, w, v, u, *z_dyn_, tol);
    ahwv.zero();
    for (size_type i=0; i<Nt_; ++i)
      ahwv.axpy(1.0, (*g_dyn_)[i]);
  }

  void applyAdjointHessian_21(Vector<Real>& ahwv,
                              const Vector<Real>& w,
                              const Vector<Real>& v,
                              const Vector<Real>& u,
                              const Vector<Real>& z,
                              Real& tol) override
  {
    for (size_type i=0; i<Nt_; ++i) {
      (*z_dyn_)[i].set(z);
      (*v_dyn_)[i].set(v);
    }
    serial_con_->applyAdjointHessian_21(ahwv, w, *v_dyn_, u, *z_dyn_, tol);
  }

  void applyAdjointHessian_22(Vector<Real>& ahwv,
                              const Vector<Real>& w,
                              const Vector<Real>& v,
                              const Vector<Real>& u,
                              const Vector<Real>& z,
                              Real& tol) override
  {
     for (size_type i=0; i<Nt_; ++i) {
      (*z_dyn_)[i].set(z);
      (*v_dyn_)[i].set(v);
    }
    serial_con_->applyAdjointHessian_22(*g_dyn_, w, *v_dyn_, u, *z_dyn_, tol);
    ahwv.zero();
    for (size_type i=0; i<Nt_; ++i)
      ahwv.axpy(1.0, (*g_dyn_)[i]);
  }
  */

  // Bring in hidden overloads
  using Constraint_SimOpt<Real>::update;
  using Constraint_SimOpt<Real>::value;
  using Constraint_SimOpt<Real>::applyAdjointJacobian_1;
  using Constraint_SimOpt<Real>::applyAdjointJacobian_2;

private:
  Ptr<SerialConstraint<Real>> serial_con_;
  size_type Nt_;
  Ptr<SerialStationaryControlsConstraintHook<Real>> hook_;
  Ptr<PartitionedVector<Real>> z_dyn_;
  Ptr<PartitionedVector<Real>> v_dyn_;
  Ptr<PartitionedVector<Real>> g_dyn_;
};
} // namespace ROL

#endif
