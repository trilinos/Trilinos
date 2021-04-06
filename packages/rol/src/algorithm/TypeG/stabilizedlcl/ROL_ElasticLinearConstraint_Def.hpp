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

#ifndef ROL_ELASTICLINEARCONSTRAINT_DEF_H
#define ROL_ELASTICLINEARCONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ElasticLinearConstraint<Real>::ElasticLinearConstraint(const Ptr<const Vector<Real>> &x,
                                                       const Ptr<Constraint<Real>>   &con,
                                                       const Ptr<const Vector<Real>> &c)
  : con_(con), x_(x->clone()), c_(c->clone()), tmp_(x->clone()) {
  setAnchor(x);
}

template<typename Real>
void ElasticLinearConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {}

template<typename Real>
void ElasticLinearConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {}

template<typename Real>
void ElasticLinearConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> xs = dynamic_cast<const PartitionedVector<Real>&>(x).get(0);
  Ptr<const Vector<Real>> xu = dynamic_cast<const PartitionedVector<Real>&>(x).get(1);
  Ptr<const Vector<Real>> xv = dynamic_cast<const PartitionedVector<Real>&>(x).get(2);
  tmp_->set(*xs); tmp_->axpy(static_cast<Real>(-1),*x_);
  con_->applyJacobian(c,*tmp_,*x_,tol);
  c.plus(*c_);
  c.plus(*xu);
  c.axpy(static_cast<Real>(-1),*xv);
}

template<typename Real>
void ElasticLinearConstraint<Real>::applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Ptr<const Vector<Real>> vs = dynamic_cast<const PartitionedVector<Real>&>(v).get(0);
  Ptr<const Vector<Real>> vu = dynamic_cast<const PartitionedVector<Real>&>(v).get(1);
  Ptr<const Vector<Real>> vv = dynamic_cast<const PartitionedVector<Real>&>(v).get(2);
  con_->applyJacobian(jv,*vs,*x_,tol);
  jv.plus(*vu);
  jv.axpy(static_cast<Real>(-1),*vv);
}

template<typename Real>
void ElasticLinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  Ptr<Vector<Real>> as = dynamic_cast<PartitionedVector<Real>&>(ajv).get(0);
  Ptr<Vector<Real>> au = dynamic_cast<PartitionedVector<Real>&>(ajv).get(1);
  Ptr<Vector<Real>> av = dynamic_cast<PartitionedVector<Real>&>(ajv).get(2);
  con_->applyAdjointJacobian(*as,v,*x_,tol);
  au->set(v.dual());
  av->set(v.dual()); av->scale(static_cast<Real>(-1));
}

template<typename Real>
void ElasticLinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol) {
  Ptr<Vector<Real>> as = dynamic_cast<PartitionedVector<Real>&>(ajv).get(0);
  Ptr<Vector<Real>> au = dynamic_cast<PartitionedVector<Real>&>(ajv).get(1);
  Ptr<Vector<Real>> av = dynamic_cast<PartitionedVector<Real>&>(ajv).get(2);
  con_->applyAdjointJacobian(*as,v,*x_,tol);
  au->set(dualv);
  av->set(dualv); av->scale(static_cast<Real>(-1));
}

template<typename Real>
void ElasticLinearConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

template<typename Real>
void ElasticLinearConstraint<Real>::setAnchor(const Ptr<const Vector<Real>> &x) {
  x_->set(*x);
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  con_->value(*c_,*x_,tol);
}

} // namespace ROL

#endif
