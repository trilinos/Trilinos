// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef L1_PENALTY_OBJECTIVE_HPP
#define L1_PENALTY_OBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "l1_utilities.hpp"

/* Reformulates the unconstrained sparse penalization problem 

   min_x f(x) + gamma ||g(x)||_1

   to the general constrained problem

   min_{x,r,s} f(x) + gamma * sum(r + s) 
   subject to  c(x,r,s) = g(x) - (r - s) = 0          
               r,s >= 0
Teuchos::GlobalMPISession::GlobalMPISession
*/

namespace ROL {

template<typename Real>
class L1PenaltyObjective : public Objective<Real> {
public: 
  L1PenaltyObjective( const Ptr<Objective<Real>>& obj,
                      const Vector<Real>&         emul,
                            Real                  gamma ) 
  : obj_(obj), e_(emul.clone()), gamma_(gamma) {
    e_->setScalar(1.0);
  }    

  void update( const Vector<Real>& x, UpdateType type, int iter = -1 ) override {
    auto [xp] = pv_cast(x);
    obj_->update(xp[X],type,iter);
  }

  Real value( const Vector<Real>& x, Real& tol ) override {
    auto [xp] = pv_cast(x);
    auto fval = obj_->value(xp[X],tol);
    return fval + gamma_*(e_->dot(xp[R])+e_->dot(xp[S]));
  }

  void gradient(       Vector<Real>& g, 
                 const Vector<Real>& x, 
                       Real&         tol ) override {
    auto [gp,xp] = pv_cast(g,x);
    obj_->gradient(gp[X],xp[X],tol);
    gp[R].setScalar(gamma_);
    gp[S].setScalar(gamma_);
  }

  void hessVec(       Vector<Real>& hv, 
                const Vector<Real>& v, 
                const Vector<Real>& x, 
                      Real&         tol ) override {
    auto [hvp,vp,xp] = pv_cast(hv,v,x);
    obj_->hessVec(hvp[X],vp[X],xp[X],tol);
    hvp[R].zero();
    hvp[S].zero();
  }

  // PartitionedVector block indices
  static constexpr int X{0}, R{1}, S{2};


private:
  Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>>    e_;
  Real gamma_;
};


} // namespace ROL



#endif //L1_PENALTY_OBJECTIVE_HPP

