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

#ifndef ROL_EXPLICIT_LINEAR_CONSTRAINT_H
#define ROL_EXPLICIT_LINEAR_CONSTRAINT_H

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_AffineTransformConstraint.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_RangeSpaceOperator.hpp"

/** @ingroup func_group
    \class ROL::ExplicitLinearConstraint
    \brief Performs null-space transformation for explicit linear equality
           constraints.

    ---
*/

namespace ROL {

template <class Real>
class ExplicitLinearConstraint {
private:
  const Ptr<Constraint<Real>> lcon_;
  const Ptr<Objective<Real>>  obj_;
  const Ptr<Vector<Real>>     x_;
  std::vector<Ptr<Constraint<Real>>> con_;

  Ptr<SimController<Real>> storage_;

  Ptr<NullSpaceOperator<Real>>        nsop_;
  Ptr<AffineTransformObjective<Real>> aobj_;
  std::vector<Ptr<AffineTransformConstraint<Real>>> acon_;

  void feasible(const Ptr<Vector<Real>> &c) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    Ptr<Vector<Real>> ran = c->clone();
    lcon_->value(*ran,*x_,tol);
    Real cnorm = ran->norm();
    if ( cnorm > static_cast<Real>(1e-4)*tol ) {
      RangeSpaceOperator<Real> rsop(lcon_,x_,c);
      Ptr<Vector<Real>> xzero = x_->clone(); xzero->zero();
      lcon_->value(*ran,*xzero,tol);
      ran->scale(static_cast<Real>(-1));
      rsop.apply(*x_,*ran,tol);
      //throw Exception::NotImplemented(">>> ExplicitLinearConstraint::feasible : Input x is not feasible!");
    }
  }

public:
  virtual ~ExplicitLinearConstraint(void) {}

  ExplicitLinearConstraint(const Ptr<Constraint<Real>> &lcon,
                           const Ptr<Objective<Real>>  &obj,
                           const Ptr<Vector<Real>>     &x,
                           const Ptr<Vector<Real>>     &c)
    : lcon_(lcon), obj_(obj), x_(x) {
    feasible(c);
    storage_ = makePtr<SimController<Real>>();
    nsop_    = makePtr<NullSpaceOperator<Real>>(lcon,x_,c);
    aobj_    = makePtr<AffineTransformObjective<Real>>(obj,nsop_,x_,storage_);
  }

  ExplicitLinearConstraint(const Ptr<Constraint<Real>> &lcon,
                           const Ptr<Objective<Real>>  &obj,
                           const Ptr<Constraint<Real>> &con,
                           const Ptr<Vector<Real>>     &x,
                           const Ptr<Vector<Real>>     &c)
    : lcon_(lcon), obj_(obj), x_(x), con_({con}) {
    feasible(c);
    storage_ = makePtr<SimController<Real>>();
    nsop_    = makePtr<NullSpaceOperator<Real>>(lcon,x_,c);
    aobj_    = makePtr<AffineTransformObjective<Real>>(obj,nsop_,x_,storage_);
    acon_.clear();
    acon_.push_back(makePtr<AffineTransformConstraint<Real>>(con,nsop_,x_,storage_));
  }

  ExplicitLinearConstraint(const Ptr<Constraint<Real>> &lcon,
                           const Ptr<Objective<Real>>  &obj,
                           const std::vector<Ptr<Constraint<Real>>> &con,
                           const Ptr<Vector<Real>>     &x,
                           const Ptr<Vector<Real>>     &c)
    : lcon_(lcon), obj_(obj), x_(x), con_(con) {
    feasible(c);
    storage_ = makePtr<SimController<Real>>();
    nsop_   = makePtr<NullSpaceOperator<Real>>(lcon,x_,c);
    aobj_   = makePtr<AffineTransformObjective<Real>>(obj,nsop_,x_,storage_);
    acon_.clear();
    int size = con_.size();
    for (int i = 0; i < size; ++i) {
      acon_.push_back(makePtr<AffineTransformConstraint<Real>>(con[i],nsop_,x_,storage_));
    }
  }

  const ROL::Ptr<Constraint<Real>> getExplicitConstraint(void) const {
    return lcon_;
  }

  const ROL::Ptr<Objective<Real>> getObjective(void) const {
    return obj_;
  }

  const ROL::Ptr<Constraint<Real>> getConstraint(const int ind = 0) const {
    if (ind < 0 || ind >= static_cast<int>(con_.size())) {
      throw Exception::NotImplemented(">>> ExplicitLinearConstraint::getConstraint : Index out of bounds!");
    }
    return con_[ind];
  }

  const ROL::Ptr<Vector<Real>> getFeasibleVector(void) const {
    return x_;
  }

  const ROL::Ptr<Objective<Real>> getTransformedObjective(void) const {
    return aobj_;
  }

  const ROL::Ptr<Constraint<Real>> getTransformedConstraint(const int ind = 0) const {
    if (ind < 0 || ind >= static_cast<int>(acon_.size())) {
      throw Exception::NotImplemented(">>> ExplicitLinearConstraint::getTransformedConstraint : Index out of bounds!");
    }
    return acon_[ind];
  }

  virtual void project(Ptr<Vector<Real>> &x,
                 const Ptr<Vector<Real>> &y) const {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    nsop_->apply(*x,*y,tol);
  }

}; // class ExplicitLinearConstraint

} // namespace ROL

#endif
