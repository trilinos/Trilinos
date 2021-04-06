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

#ifndef ROL_CAUCHYPOINT_U_H
#define ROL_CAUCHYPOINT_U_H

/** \class ROL::CauchyPoint_U
    \brief Provides interface for the Cauchy point trust-region subproblem solver.
*/

#include "ROL_TrustRegion_U.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class CauchyPoint_U : public TrustRegion_U<Real> {
private:

  Ptr<Vector<Real>> dual_;

public:

  CauchyPoint_U() {}

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    dual_ = g.clone();
  }

  void solve( Vector<Real> &s, Real &snorm, Real &pRed,
              int &iflag, int &iter, const Real del,
              TrustRegionModel_U<Real> &model) {
    const Real zero(0), half(0.5);
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Set step to (projected) gradient
    s.set(model.getGradient()->dual());
    // Apply (reduced) Hessian to (projected) gradient
    model.hessVec(*dual_,s,s,tol);
    Real gnorm  = s.norm();
    Real gnorm2 = gnorm*gnorm;
    //Real gBg    = dual_->dot(s.dual());
    Real gBg    = dual_->apply(s);
    Real alpha  = gnorm2/gBg;
    if ( alpha*gnorm >= del || gBg <= zero ) {
      alpha = del/gnorm;
    }
    s.scale(-alpha);
    snorm = alpha*gnorm;
    iflag = 0;
    iter  = 0;
    pRed  = alpha*(gnorm2 - half*alpha*gBg);
  }
};

} // namespace ROL

#endif
