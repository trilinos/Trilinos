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

#ifndef ROL_TRUNCATEDCG_U_H
#define ROL_TRUNCATEDCG_U_H

/** \class ROL::TruncatedCG_U
    \brief Provides interface for truncated CG trust-region subproblem solver.
*/

#include "ROL_TrustRegion_U.hpp"
#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class TruncatedCG_U : public TrustRegion_U<Real> {
private:
  Ptr<Vector<Real>> s_, g_, v_, p_, Hp_;

  int maxit_;
  Real tol1_;
  Real tol2_;

public:

  // Constructor
  TruncatedCG_U( ParameterList &parlist ) {
    // Unravel Parameter List
    Real em4(1e-4), em2(1e-2);
    maxit_ = parlist.sublist("General").sublist("Krylov").get("Iteration Limit",20);
    tol1_  = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance",em4);
    tol2_  = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance",em2);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    s_ = x.clone();
    v_ = x.clone();
    p_ = x.clone();
    g_ = g.clone();
    Hp_ = g.clone();
  }

  void solve( Vector<Real>             &s,
              Real                     &snorm,
              Real                     &pRed,
              int                      &iflag,
              int                      &iter,
              const Real                del,
              TrustRegionModel_U<Real> &model ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    const Real zero(0), one(1), two(2), half(0.5);
    // Initialize step
    s.zero(); s_->zero();
    snorm = zero;
    Real snorm2(0), s1norm2(0);
    // Compute (projected) gradient
    g_->set(*model.getGradient());
    Real gnorm = g_->norm(), normg = gnorm;
    const Real gtol = std::min(tol1_,tol2_*gnorm);
    // Preconditioned (projected) gradient vector
    model.precond(*v_,*g_,s,tol);
    // Initialize basis vector
    p_->set(*v_); p_->scale(-one);
    //Real pnorm2 = v_->dot(g_->dual());
    Real pnorm2 = v_->apply(*g_);
    if ( pnorm2 <= zero ) {
      iflag = 4;
      iter  = 0;
      return;
    }
    // Initialize scalar storage
    iter = 0; iflag = 0;
    Real kappa(0), beta(0), sigma(0), alpha(0), tmp(0), sMp(0);
    Real gv = pnorm2; //v_->dot(g_->dual());
    pRed = zero;
    // Iterate CG
    for (iter = 0; iter < maxit_; iter++) {
      // Apply Hessian to direction p
      model.hessVec(*Hp_,*p_,s,tol);
      // Check positivity of Hessian
      //kappa = p_->dot(Hp_->dual());
      kappa = p_->apply(*Hp_);
      if (kappa <= zero) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p_);
        snorm = del;
        iflag = 2;
        break;
      }
      // Update step
      alpha = gv/kappa;
      s_->set(s);
      s_->axpy(alpha,*p_);
      s1norm2 = snorm2 + two*alpha*sMp + alpha*alpha*pnorm2;
      // Check if step exceeds trust region radius
      if (s1norm2 >= del*del) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p_);
        snorm = del;
        iflag = 3;
        break;
      }
      // Update model predicted reduction
      pRed += half*alpha*gv;
      // Set step to temporary step and store norm
      s.set(*s_);
      snorm2 = s1norm2;
      // Check for convergence
      g_->axpy(alpha,*Hp_);
      normg = g_->norm();
      if (normg < gtol) {
        break;
      }
      // Preconditioned updated (projected) gradient vector
      model.precond(*v_,*g_,s,tol);
      tmp   = gv;
      //gv    = v_->dot(g_->dual());
      gv    = v_->apply(*g_);
      beta  = gv/tmp;
      // Update basis vector
      p_->scale(beta);
      p_->axpy(-one,*v_);
      sMp    = beta*(sMp+alpha*pnorm2);
      pnorm2 = gv + beta*beta*pnorm2;
    }
    // Update model predicted reduction
    if (iflag > 0) {
      pRed += sigma*(gv-half*sigma*kappa);
    }
    else {
      snorm = std::sqrt(snorm2);
    }
    // Check iteration count
    if (iter == maxit_) {
      iflag = 1;
    }
    if (iflag != 1) {
      iter++;
    }
  }
};

} // namespace ROL

#endif
