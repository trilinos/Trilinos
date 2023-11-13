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

