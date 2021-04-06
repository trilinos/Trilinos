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

#ifndef ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_DEF_H
#define ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ScalarLinearConstraint<Real>::ScalarLinearConstraint(const Ptr<const Vector<Real>> &a,
                                                     const Real b)
  : a_(a), b_(b) {}

template<typename Real>
void ScalarLinearConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  SingletonVector<Real> &cc = dynamic_cast<SingletonVector<Real>&>(c);
  //cc.setValue(a_->dot(x.dual()) - b_);
  cc.setValue(a_->apply(x) - b_);
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                           const Vector<Real> &v,
                                           const Vector<Real> &x, Real &tol) {
  SingletonVector<Real> &jc = dynamic_cast<SingletonVector<Real>&>(jv);
  //jc.setValue(a_->dot(v.dual()));
  jc.setValue(a_->apply(v));
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                                  const Vector<Real> &v,
                                                  const Vector<Real> &x, Real &tol) {
  const SingletonVector<Real>&    vc = dynamic_cast<const SingletonVector<Real>&>(v);
  ajv.set(*a_);
  ajv.scale(vc.getValue());
}

template<typename Real>
void ScalarLinearConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &x, Real &tol) {
  ahuv.zero();
}

template<typename Real>
std::vector<Real> ScalarLinearConstraint<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                                     Vector<Real> &v2,
                                                               const Vector<Real> &b1,
                                                               const Vector<Real> &b2,
                                                               const Vector<Real> &x, Real &tol) {
  SingletonVector<Real>&       v2c = dynamic_cast<SingletonVector<Real>&>(v2);
  const SingletonVector<Real>& b2c = dynamic_cast<const SingletonVector<Real>&>(b2);

  //v2c.setValue( (a_->dot(b1.dual()) - b2c.getValue() )/a_->dot(*a_) );
  v2c.setValue( (a_->apply(b1) - b2c.getValue() )/a_->dot(*a_) );
  v1.set(b1.dual());
  v1.axpy(-v2c.getValue(),a_->dual());

  std::vector<Real> out;
  return out;
}

} // namespace ROL

#endif
