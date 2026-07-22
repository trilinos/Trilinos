// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTED_SYSTEM_OPERATOR_H
#define ROL_AUGMENTED_SYSTEM_OPERATOR_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::AugmentedSystemOperator
    \brief Apply the augmented system operator.

    ---
*/

namespace ROL {

template <class Real>
class AugmentedSystemOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>>   con_;
  const Ptr<const Vector<Real>> x_;
  const Real                    delta_;

public:
  virtual ~AugmentedSystemOperator() {}
  AugmentedSystemOperator(const Ptr<Constraint<Real>>   &con,
                          const Ptr<const Vector<Real>> &x,
                          const Real                     delta = 0.0)
    : con_(con), x_(x), delta_(delta) {}

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    PartitionedVector<Real>      &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
    const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

    Ptr<Vector<Real>> h1 = Hvp.get(0)->dual().clone();
    con_->applyAdjointJacobian(*h1, *(vp.get(1)), *x_, tol);
    //con_->applyAdjointJacobian(*(Hvp.get(0)), *(vp.get(1)), *x_, tol);
    //Hvp.get(0)->plus(*(vp.get(0)));
    Hvp.get(0)->set(h1->dual()); Hvp.get(0)->plus(*(vp.get(0)));

    con_->applyJacobian(*(Hvp.get(1)), *(vp.get(0)), *x_, tol);
    if ( delta_ > static_cast<Real>(0) ) {
      Hvp.get(1)->axpy(-delta_*delta_, *(vp.get(1)));
    }
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    apply(Hv,v,tol);
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemOperator::applyInverse : Not implemented!");
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemOperator::applyAdjointInverse : Not implemented!");
  }

}; // class AugmentedSystemOperator

} // namespace ROL

#endif
