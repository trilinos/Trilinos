// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_L1OBJECTIVE_H
#define ROL_L1OBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::l1Objective
    \brief Provides the interface to evaluate the weighted/shifted l1 objective function.

        ---
*/


namespace ROL {

template<typename Real>
class l1Objective : public Objective<Real> {
private:
  const Ptr<Vector<Real>> weights_, shift_;
  Ptr<Vector<Real>> tmp_;

  struct ProjSymBnd : public Elementwise::BinaryFunction<Real> {
       Real apply(const Real &xc, const Real &yc) const { return std::min(yc, std::max(-yc, xc)); }
  } psb_;
 
public:

  l1Objective(const Ptr<Vector<Real>> &weights)
    : weights_(weights), shift_(weights->dual().clone()) {
    shift_->zero();
    tmp_ = shift_->clone();
  }

  l1Objective(const Ptr<Vector<Real>> &weights, const Ptr<Vector<Real>> &shift)
    : weights_(weights),  shift_(shift) {
    tmp_ = shift_->clone();
  }
  
  Real value( const Vector<Real> &x, Real &tol ) {
    tmp_->set(x);
    tmp_->axpy(static_cast<Real>(-1),*shift_);
    tmp_->applyUnary(Elementwise::AbsoluteValue<Real>());
    return weights_->apply(*tmp_);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.set(x);
    g.axpy(static_cast<Real>(-1),*shift_);
    g.applyUnary(Elementwise::Sign<Real>());
    g.applyBinary(Elementwise::Multiply<Real>(), *weights_);
  }

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    gradient(*tmp_, x, tol);
    return tmp_->apply(d);
  }

  void prox( Vector<Real> &Pv, const Vector<Real> &v, Real t, Real &tol){
    Pv.set(*shift_);
    Pv.axpy(static_cast<Real>(-1), v);
    Pv.scale(static_cast<Real>(1) / t);
    Pv.applyBinary(psb_, *weights_);
    Pv.scale(t);
    Pv.plus(v);
  } 
}; // class l1Objective

} // namespace ROL

#endif
