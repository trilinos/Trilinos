// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef L1_PENALTY_CONSTRAINT_HPP
#define L1_PENALTY_CONSTRAINT_HPP

#include "ROL_Constraint.hpp"
#include "l1_utilities.hpp"

namespace ROL {

template<typename Real>
class L1PenaltyConstraint : public Constraint<Real> {
public:

  L1PenaltyConstraint( const Ptr<Constraint<Real>>& con ) 
  : con_(con) {}

  void update( const Vector<Real>& x, 
                     UpdateType    type, 
                     int           iter = -1 ) override {
    auto [xp] = pv_cast(x);
    con_->update(xp[X],type,iter);
  }

  void value(       Vector<Real>& c,
              const Vector<Real>& x,
                    Real&         tol ) override {
    auto [xp] = pv_cast(x);
    con_->value(c,xp[X],tol);
    c.plus(xp[S]);
    c.axpy(-1.0,xp[R]);
  }

  void applyJacobian(       Vector<Real>& jv,
                      const Vector<Real>& v,
                      const Vector<Real>& x,
                            Real&         tol ) override {
    auto [vp,xp] = pv_cast(v,x);
    con_->applyJacobian(jv,vp[X],xp[X],tol);
    jv.plus(vp[S]);
    jv.axpy(-1.0,vp[R]);
  }

  void applyAdjointJacobian(       Vector<Real>& ajv,
                             const Vector<Real>& v,
                             const Vector<Real>& x,
                                   Real&         tol ) override {
    auto [ajvp,xp] = pv_cast(ajv,x);  
    con_->applyAdjointJacobian(ajvp[X],v,xp[X],tol);
    ajvp[R].set(v);
    ajvp[R].scale(-1);
    ajvp[S].set(v);
  }

  void applyAdjointHessian(      Vector<Real>& ahuv,
                           const Vector<Real>& u,
                           const Vector<Real>& v,
                           const Vector<Real>& x,
                                 Real&         tol ) override {
    auto [ahuvp,vp,xp] = pv_cast(ahuv,v,x);
    con_->applyAdjointHessian(ahuvp[X],u,vp[X],xp[X],tol);
    ahuvp[R].zero();
    ahuvp[S].zero();
  }

  // PartitionedVector block indices
  static constexpr int X{0}, R{1}, S{2};

private:

  Ptr<Constraint<Real>> con_;

};


} // namespace ROL

#endif //L1_PENALTY_CONSTRAINT_HPP

