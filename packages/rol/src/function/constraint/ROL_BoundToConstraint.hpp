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

#ifndef ROL_BOUND_TO_CONSTRAINT_H
#define ROL_BOUND_TO_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_LowerBoundToConstraint.hpp"
#include "ROL_UpperBoundToConstraint.hpp"

/**  @ingroup func_group
     \class ROL::BoundToConstraint 
     \brief Provides an implementation for bound constraints.
*/

namespace ROL {

template <class Real>
class BoundToConstraint : public Constraint<Real> { 
private:
  ROL::Ptr<InequalityConstraint<Real> > lo_;
  ROL::Ptr<InequalityConstraint<Real> > up_;
  ROL::Ptr<Vector<Real> > tmp_;

public:
  BoundToConstraint(BoundConstraint<Real> &bnd, const Vector<Real> &x) {
    lo_ = ROL::makePtr<LowerBoundToConstraint<Real>>(bnd,x);
    up_ = ROL::makePtr<UpperBoundToConstraint<Real>>(bnd,x);
    tmp_ = x.clone();
  }

  BoundToConstraint(const Vector<Real> &lo, const Vector<Real> &up) {
    lo_ = ROL::makePtr<LowerBoundToConstraint<Real>>(lo);
    up_ = ROL::makePtr<UpperBoundToConstraint<Real>>(up);
    tmp_ = lo.clone();
  }

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    Vector<Real> &c0 = *(dynamic_cast<PartitionedVector<Real>&>(c).get(0));
    Vector<Real> &c1 = *(dynamic_cast<PartitionedVector<Real>&>(c).get(1));
    lo_->value(c0,x,tol);
    up_->value(c1,x,tol);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    Vector<Real> &jv0 = *(dynamic_cast<PartitionedVector<Real>&>(jv).get(0));
    Vector<Real> &jv1 = *(dynamic_cast<PartitionedVector<Real>&>(jv).get(1));
    lo_->applyJacobian(jv0,v,x,tol);
    up_->applyJacobian(jv1,v,x,tol);
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    const Vector<Real> &v0 = *(dynamic_cast<const PartitionedVector<Real>&>(v).get(0));
    const Vector<Real> &v1 = *(dynamic_cast<const PartitionedVector<Real>&>(v).get(1));
    lo_->applyAdjointJacobian(ajv,v0,x,tol);
    up_->applyAdjointJacobian(*tmp_,v1,x,tol);
    ajv.plus(*tmp_); 
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {
    ahuv.zero();
  }
};

}

#endif
