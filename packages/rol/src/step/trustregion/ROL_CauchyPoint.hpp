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

#ifndef ROL_CAUCHYPOINT_H
#define ROL_CAUCHYPOINT_H

/** \class ROL::CauchyPoint
    \brief Provides interface for the Cauchy point trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL { 

template<class Real>
class CauchyPoint : public TrustRegion<Real> {
private:

  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > p_;
  Teuchos::RCP<Vector<Real> > Hp_;

  Real pRed_;
  Real eps_;
  Real alpha_;

  bool useCGTCP_;

public:

  // Constructor
  CauchyPoint( Teuchos::ParameterList &parlist ) 
    : TrustRegion<Real>(parlist), pRed_(0.0), alpha_(-1.0), useCGTCP_(false) {
    // Unravel Parameter List
    Real TRsafe = parlist.get("Trust-Region Safeguard",100.0);
    eps_        = TRsafe*ROL_EPSILON;
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    Hp_ = g.clone();
    if ( useCGTCP_ ) {
      g_ = g.clone();
      p_ = s.clone();
    }
  }

  void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
            const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) { 
    if ( pObj.isConActivated() ) {
      if ( useCGTCP_ ) {
        cauchypoint_CGT( s, snorm, del, iflag, iter, x, grad, gnorm, pObj );
      } 
      else {
        cauchypoint_M( s, snorm, del, iflag, iter, x, grad, gnorm, pObj );
      }
    }
    else {
      cauchypoint_unc( s, snorm, del, iflag, iter, x, grad, gnorm, pObj );
    }
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }

private:
  void cauchypoint_unc( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                        const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);
    pObj.hessVec(*Hp_,grad.dual(),x,tol);
    Real gBg = Hp_->dot(grad);
    Real tau = 1.0;
    if ( gBg > 0.0 ) {
      tau = std::min(1.0, gnorm*gnorm*gnorm/gBg);
    }

    s.set(grad.dual());
    s.scale(-tau*del/gnorm);
    snorm = tau*del;
    iflag = 0;
    iter  = 0;
    pRed_ = tau*del/gnorm * pow(gnorm,2.0) - 0.5*pow(tau*del/gnorm,2.0)*gBg;
  }

  void cauchypoint_M( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                      const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Parameters
    Real mu0   = 1.e-2;
    Real mu1   = 1.0;
    Real beta1 = 0.0;
    Real beta2 = 0.0;
    bool decr  = true;
    bool stat  = true;

    // Initial step length
    Real alpha  = 1.0;
    if ( alpha_ > 0.0 ) {
      alpha = alpha_;
    } 
    Real alpha0   = alpha;
    Real alphamax = 1.e4*alpha;
    
    // Initial model value
    s.set(grad.dual());
    s.scale(-alpha);
    pObj.computeProjectedStep(s,x);
    snorm = s.norm();
    pObj.hessVec(*Hp_,s,x,tol);
    Real gs   = s.dot(grad.dual());
    Real val  = gs + 0.5*s.dot(Hp_->dual());
    Real val0 = val;

    // Determine whether to increase or decrease alpha
    if ( val > mu0 * gs || snorm > mu1 * del ) { 
      beta1 = 0.5; 
      beta2 = 0.5; 
      decr  = true;
    }
    else {
      beta1 = 2.0;
      beta2 = 2.0;
      decr  = false;
    }

    while ( stat ) {
      // Update step length
      alpha0 = alpha;
      val0   = val;
      alpha *= (beta1+beta2)*0.5;
  
      // Update model value
      s.set(grad.dual());
      s.scale(-alpha);
      pObj.computeProjectedStep(s,x);
      snorm = s.norm();
      pObj.hessVec(*Hp_,s,x,tol);
      gs    = s.dot(grad.dual());
      val   = gs + 0.5*s.dot(Hp_->dual());

      // Update termination criterion
      if ( decr ) {
        stat = ( val > mu0 * gs || snorm > mu1 * del );
        if ( std::abs(val) < eps_ && std::abs(mu0 *gs) < eps_ ) {
          stat = (snorm > mu1 * del);
        }
      }
      else {
        stat = !( val > mu0 * gs || snorm > mu1 * del );
        if ( std::abs(val) < eps_ && std::abs(mu0 *gs) < eps_ ) {
          stat = !(snorm > mu1 * del);
        }
        if ( alpha > alphamax ) {
          stat = false;
        }
      } 
    }
    // Reset to last 'successful' step
    val   = val0;
    alpha = alpha0;
    s.set(grad.dual());
    s.scale(-alpha);
    pObj.computeProjectedStep(s,x);
    snorm = s.norm();
    
    alpha_ = alpha;
    pRed_  = -val;
  }

  void cauchypoint_CGT( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                        const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);
    bool tmax_flag = true;
    int maxit      = 20;
    Real t         = del/gnorm;
    Real tmax      = 1.e10;
    Real tmin      = 0.0;
    Real gs        = 0.0;
    Real c1        = 0.25;
    Real c2        = 0.75;
    Real c3        = 0.9;
    Real c4        = 0.25;
    Real pgnorm    = 0.0;
    for ( int i = 0; i < maxit; i++ ) {
      // Compute p = x + s = P(x - t*g)
      p_->set(x);
      p_->axpy(-t,grad.dual()); 
      pObj.project(*p_);
      // Compute s = p - x = P(x - t*g) - x
      s.set(*p_);
      s.axpy(-1.0,x);
      snorm = s.norm();
      // Evaluate Model
      pObj.hessVec(*Hp_,s,x,tol);
      gs = s.dot(grad.dual());
      pRed_ = -gs - 0.5*s.dot(Hp_->dual());

      // Check Stopping Conditions
      g_->set(grad);
      pObj.pruneActive(*g_,grad,*p_); // Project gradient onto tangent cone at p
      pgnorm = g_->norm();
      if ( snorm > del || pRed_ < -c2*gs ) {
        tmax = t;
        tmax_flag = false;
      }
      else if ( snorm < c3*del && pRed_ > -c1*gs && pgnorm > c4*std::abs(gs)/del ) {
        tmin = t;
      } 
      else {
        break;
      }
   
      // Update t
      if ( tmax_flag ) {
        t *= 2.0;
      }
      else {
        t = 0.5*(tmax + tmin);
      }
    }
  }
};

}

#endif
