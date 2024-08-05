// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
