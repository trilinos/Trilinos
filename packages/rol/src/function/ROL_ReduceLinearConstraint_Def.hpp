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

#ifndef ROL_REDUCE_LINEAR_CONSTRAINT_DEF_H
#define ROL_REDUCE_LINEAR_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
ReduceLinearConstraint<Real>::ReduceLinearConstraint(const Ptr<Constraint<Real>>   &lcon,
                                                     const Ptr<Vector<Real>>       &x,
                                                     const Ptr<const Vector<Real>> &c)
  : lcon_(lcon), x_(x),
    storage_(makePtr<VectorController<Real>>()),
    nsop_(makePtr<NullSpaceOperator<Real>>(lcon,x_,c)) {
  feasible(c);
}

template<typename Real>
Ptr<Objective<Real>> ReduceLinearConstraint<Real>::transform(const Ptr<Objective<Real>> &obj) const {
  return makePtr<AffineTransformObjective<Real>>(obj,nsop_,x_,storage_);
}

template<typename Real>
Ptr<Constraint<Real>> ReduceLinearConstraint<Real>::transform(const Ptr<Constraint<Real>> &con) const {
  return makePtr<AffineTransformConstraint<Real>>(con,nsop_,x_,storage_);
}

template<typename Real>
Ptr<Constraint<Real>> ReduceLinearConstraint<Real>::getLinearConstraint(void) const {
  return lcon_;
}

template<typename Real>
Ptr<const Vector<Real>> ReduceLinearConstraint<Real>::getFeasibleVector(void) const {
  return x_;
}

template<typename Real>
void ReduceLinearConstraint<Real>::project(Vector<Real> &x, const Vector<Real> &y) const {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  nsop_->apply(x,y,tol);
}

template<typename Real>
void ReduceLinearConstraint<Real>::project(const Ptr<Vector<Real>> &x, const Ptr<const Vector<Real>> &y) const {
  project(*x,*y);
}

template<typename Real>
void ReduceLinearConstraint<Real>::feasible(const Ptr<const Vector<Real>> &c) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  Ptr<Vector<Real>> ran = c->clone();
  lcon_->value(*ran,*x_,tol);
  Real cnorm = ran->norm();
  if ( cnorm > static_cast<Real>(1e-4)*tol ) {
    RangeSpaceOperator<Real> rsop(lcon_,x_,c);
    Ptr<Vector<Real>> xzero = x_->clone(); xzero->zero();
    lcon_->value(*ran,*xzero,tol);
    ran->scale(static_cast<Real>(-1));
    nsop_->apply(*xzero,*x_,tol);
    rsop.apply(*x_,*ran,tol);
    x_->plus(*xzero);
    //throw Exception::NotImplemented(">>> ReduceLinearConstraint::feasible : Input x is not feasible!");
  }
}

} // namespace ROL

#endif
