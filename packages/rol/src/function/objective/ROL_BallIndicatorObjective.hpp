// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BALLINDICATOROBJECTIVE_H
#define ROL_BALLINDICATOROBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::BallIndicatorObjective
    \brief Provides the interface to evaluate the indicator function of norm constraints.

        ---
*/


namespace ROL {

template<typename Real>
class BallIndicatorObjective : public Objective<Real> {
private:
  const Ptr<Vector<Real>> x_, pwa_;
  const Real rad_;
 
public:

  BallIndicatorObjective(const Ptr<Vector<Real>> &x, Real rad)
    : x_(x), pwa_(x->clone()), rad_(rad) {}

  Real value( const Vector<Real> &x, Real &tol ) {
    const Real zero(0), one(1);
    pwa_->set(x); pwa_->axpy(-one,*x_);
    Real norm = pwa_->norm();
    return (norm <= rad_) ? zero : ROL_INF<Real>();
  }

  void prox( Vector<Real> &Pv, const Vector<Real> &v, Real t, Real &tol){
    pwa_->set(v); pwa_->axpy(-one,*x_);
    Real norm = pwa_->norm();
    if(norm <= rad_) {
      Pv.set(v);
    }
    else {
      Pv.set(*x_);
      Pv.axpy(rad_/norm,*pwa_);
    }
  }
}; // class BallIndicatorObjective

} // namespace ROL

#endif
