// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AUGMENTED_SYSTEM_PREC_OPERATOR_H
#define ROL_AUGMENTED_SYSTEM_PREC_OPERATOR_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::AugmentedSystemPrecOperator
    \brief Implements a preconditioner for the augmented system.

    ---
*/

namespace ROL {

template <class Real>
class AugmentedSystemPrecOperator : public LinearOperator<Real> {
private:
  const Ptr<Constraint<Real>>   con_;
  const Ptr<const Vector<Real>> x_;

public:
  virtual ~AugmentedSystemPrecOperator() {}
  AugmentedSystemPrecOperator(const Ptr<Constraint<Real>>   &con,
                              const Ptr<const Vector<Real>> &x)
    : con_(con), x_(x) {}

  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemPrecOperator::apply : Not Implemented!");
  }

  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    throw Exception::NotImplemented(">>> AugmentedSystemPrecOperator::applyAdjoint : Not Implemented!");
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Real zero(0);
    PartitionedVector<Real>      &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
    const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);

    Hvp.set(0, *(vp.get(0)));
    // Second x should be dual, but unused?
    con_->applyPreconditioner(*(Hvp.get(1)),*(vp.get(1)),*x_,*x_, zero);
  }

  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    applyInverse(Hv,v,tol);
  }

}; // class AugmentedSystemPrecOperator

} // namespace ROL

#endif
