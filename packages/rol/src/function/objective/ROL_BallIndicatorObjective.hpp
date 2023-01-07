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
