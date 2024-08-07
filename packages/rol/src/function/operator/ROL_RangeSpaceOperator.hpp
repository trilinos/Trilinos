// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RANGE_SPACE_OPERATOR_H
#define ROL_RANGE_SPACE_OPERATOR_H

#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::RangeSpaceOperator
    \brief Projects on to the null space of a linear constraint.

    ---
*/

namespace ROL {

template <class Real>
class RangeSpaceOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>> con_;
  const Ptr<Vector<Real>> x_;

  mutable Ptr<Vector<Real>> b1_;
  mutable Ptr<Vector<Real>> b2_;
  mutable Ptr<Vector<Real>> mul_;

public:
  virtual ~RangeSpaceOperator() {}
  RangeSpaceOperator(const Ptr<Constraint<Real>>   &con,
                     const Ptr<const Vector<Real>> &dom,
                     const Ptr<const Vector<Real>> &ran)
    : con_(con), x_(dom->clone()) {
    x_->set(*dom);
    b1_  = dom->dual().clone();
    b2_  = ran->clone();
    mul_ = ran->dual().clone();
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    x_->set(x);
  }

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    b1_->zero(); b2_->set(v);
    con_->solveAugmentedSystem(Hv,*mul_,*b1_,*b2_,*x_,tol); // This assumes linearity
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> RangeSpaceOperator::applyAdjoint : Not Implemented!");
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> RangeSpaceOperator::applyInverse : Not Implemented!");
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> RangeSpaceOperator::applyAdjointInverse : Not Implemented!");
  }

}; // class RangeSpaceOperator

} // namespace ROL

#endif
