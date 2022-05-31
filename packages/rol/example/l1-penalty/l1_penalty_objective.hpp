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

