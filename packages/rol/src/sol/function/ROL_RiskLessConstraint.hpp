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

#ifndef ROL_RISKLESS_CONSTRAINT_H
#define ROL_RISKLESS_CONSTRAINT_H

#include "ROL_RiskVector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

template <class Real>
class RiskLessConstraint : public Constraint<Real> {
private:
  const ROL::Ptr<Constraint<Real> > con_;

public:
  RiskLessConstraint(const ROL::Ptr<Constraint<Real> > &con)
    : con_(con) {}

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    ROL::Ptr<const Vector<Real> > x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    con_->value(c,*x0,tol);
  }

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    ROL::Ptr<const Vector<Real> > x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    ROL::Ptr<const Vector<Real> > v0 = dynamic_cast<const RiskVector<Real>&>(v).getVector();
    con_->applyJacobian(jv,*v0,*x0,tol);
  }

  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    ROL::Ptr<const Vector<Real> > x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    ROL::Ptr<Vector<Real> > ajv0 = dynamic_cast<RiskVector<Real>&>(ajv).getVector();
    con_->applyAdjointJacobian(*ajv0,v,*x0,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    ROL::Ptr<const Vector<Real> > x0 = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    ROL::Ptr<const Vector<Real> > v0 = dynamic_cast<const RiskVector<Real>&>(v).getVector();
    ROL::Ptr<Vector<Real> > ahuv0 = dynamic_cast<RiskVector<Real>&>(ahuv).getVector();
    con_->applyAdjointHessian(*ahuv0,u,*v0,*x0,tol);
  }

public:
  virtual void setParameter(const std::vector<Real> &param) {
    Constraint<Real>::setParameter(param);
    con_->setParameter(param);
  }

}; // class RiskLessConstraint

} // namespace ROL

#endif
