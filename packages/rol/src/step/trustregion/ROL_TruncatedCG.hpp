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

#ifndef ROL_TRUNCATEDCG_H
#define ROL_TRUNCATEDCG_H

/** \class ROL::TruncatedCG
    \brief Provides interface for truncated CG trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

namespace ROL { 

template<class Real>
class TruncatedCG : public TrustRegion<Real> {
private:

  Teuchos::RCP<Vector<Real> > s_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > v_;
  Teuchos::RCP<Vector<Real> > p_;
  Teuchos::RCP<Vector<Real> > Hp_;

  int maxit_;
  Real tol1_;
  Real tol2_;

  Real pRed_;

public:

  // Constructor
  TruncatedCG( Teuchos::ParameterList &parlist ) : TrustRegion<Real>(parlist), pRed_(0.0) {
    // Unravel Parameter List
    maxit_ = parlist.get("Maximum Number of Krylov Iterations", 20);
    tol1_  = parlist.get("Absolute Krylov Tolerance",           1.e-4);
    tol2_  = parlist.get("Relative Krylov Tolerance",           1.e-2);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    s_ = s.clone();
    g_ = g.clone();
    v_ = s.clone();
    p_ = s.clone();
    Hp_ = g.clone();
  }

  void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
            const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) { 
    Real tol = std::sqrt(ROL_EPSILON);
    const Real gtol = std::min(tol1_,tol2_*gnorm);

    // Old and New Step Vectors
    s.zero(); 
    snorm = 0.0;
    Real snorm2  = 0.0;
    s_->zero();
    Real s1norm2 = 0.0;

    // Gradient Vector
    g_->set(grad);
    Real normg = gnorm;
    if ( pObj.isConActivated() ) {
      pObj.pruneActive(*g_,grad,x);
      normg = g_->norm();
    }

    // Preconditioned Gradient Vector
    //pObj.precond(*v,*g,x,tol);
    pObj.reducedPrecond(*v_,*g_,x,grad,x,tol);

    // Basis Vector
    p_->set(*v_); 
    p_->scale(-1.0);
    Real pnorm2 = v_->dot(g_->dual());

    iter        = 0; 
    iflag       = 0;
    Real kappa  = 0.0;
    Real beta   = 0.0; 
    Real sigma  = 0.0; 
    Real alpha  = 0.0; 
    Real tmp    = 0.0;
    Real gv     = v_->dot(g_->dual());
    Real sMp    = 0.0;
    pRed_       = 0.0;

    for (iter = 0; iter < maxit_; iter++) {
      //pObj.hessVec(*Hp,*p,x,tol);
      pObj.reducedHessVec(*Hp_,*p_,x,grad,x,tol);

      kappa = p_->dot(Hp_->dual());
      if (kappa <= 0.0) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p_);
        iflag = 2; 
        break;
      }

      alpha = gv/kappa;
      s_->set(s); 
      s_->axpy(alpha,*p_);
      s1norm2 = snorm2 + 2.0*alpha*sMp + alpha*alpha*pnorm2;

      if (s1norm2 >= del*del) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p_);
        iflag = 3; 
        break;
      }

      pRed_ += 0.5*alpha*gv;

      s.set(*s_);
      snorm2 = s1norm2;  

      g_->axpy(alpha,*Hp_);
      normg = g_->norm();
      if (normg < gtol) {
        break;
      }

      //pObj.precond(*v,*g,x,tol);
      pObj.reducedPrecond(*v_,*g_,x,grad,x,tol);
      tmp   = gv; 
      gv    = v_->dot(g_->dual());
      beta  = gv/tmp;    

      p_->scale(beta);
      p_->axpy(-1.0,*v_);
      sMp    = beta*(sMp+alpha*pnorm2);
      pnorm2 = gv + beta*beta*pnorm2; 
    }
    if (iflag > 0) {
      pRed_ += sigma*(gv-0.5*sigma*kappa);
    }

    if (iter == maxit_) {
      iflag = 1;
    }
    if (iflag != 1) { 
      iter++;
    }

    snorm = s.norm();
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }

#if 0
  void truncatedCG_proj( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                         const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);

    const Real gtol = std::min(tol1_,tol2_*gnorm);

    // Compute Cauchy Point
    Real scnorm = 0.0;
    Teuchos::RCP<Vector<Real> > sc = x.clone();
    cauchypoint(*sc,scnorm,del,iflag,iter,x,grad,gnorm,pObj);
    Teuchos::RCP<Vector<Real> > xc = x.clone();
    xc->set(x);
    xc->plus(*sc);

    // Old and New Step Vectors
    s.set(*sc); 
    snorm = s.norm();
    Real snorm2  = snorm*snorm;
    Teuchos::RCP<Vector<Real> > s1 = x.clone();
    s1->zero();
    Real s1norm2 = 0.0;

    // Gradient Vector
    Teuchos::RCP<Vector<Real> > g = x.clone(); 
    g->set(grad);
    Teuchos::RCP<Vector<Real> > Hs = x.clone();
    pObj.reducedHessVec(*Hs,s,*xc,x,tol);
    g->plus(*Hs);
    Real normg = g->norm();

    // Preconditioned Gradient Vector
    Teuchos::RCP<Vector<Real> > v  = x.clone();
    pObj.reducedPrecond(*v,*g,*xc,x,tol);

    // Basis Vector
    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 
    p->scale(-1.0);
    Real pnorm2 = v->dot(*g);

    // Hessian Times Basis Vector
    Teuchos::RCP<Vector<Real> > Hp = x.clone();

    iter        = 0; 
    iflag       = 0; 
    Real kappa  = 0.0;
    Real beta   = 0.0; 
    Real sigma  = 0.0; 
    Real alpha  = 0.0; 
    Real tmp    = 0.0;
    Real gv     = v->dot(*g);
    Real sMp    = 0.0;
    pRed_ = 0.0;

    for (iter = 0; iter < maxit_; iter++) {
      pObj.reducedHessVec(*Hp,*p,*xc,x,tol);

      kappa = p->dot(*Hp);
      if (kappa <= 0) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p);
        iflag = 2; 
        break;
      }

      alpha = gv/kappa;
      s1->set(s); 
      s1->axpy(alpha,*p);
      s1norm2 = snorm2 + 2.0*alpha*sMp + alpha*alpha*pnorm2;

      if (s1norm2 >= del*del) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p);
        iflag = 3; 
        break;
      }

      pRed_ += 0.5*alpha*gv;

      s.set(*s1);
      snorm2 = s1norm2;  

      g->axpy(alpha,*Hp);
      normg = g->norm();
      if (normg < gtol) {
        break;
      }

      pObj.reducedPrecond(*v,*g,*xc,x,tol);
      tmp   = gv; 
      gv    = v->dot(*g);
      beta  = gv/tmp;    

      p->scale(beta);
      p->axpy(-1.0,*v);
      sMp    = beta*(sMp+alpha*pnorm2);
      pnorm2 = gv + beta*beta*pnorm2; 
    }
    if (iflag > 0) {
      pRed_ += sigma*(gv-0.5*sigma*kappa);
    }

    if (iter == maxit_) {
      iflag = 1;
    }
    if (iflag != 1) { 
      iter++;
    }

    snorm = s.norm();
  }
#endif
};

}

#endif
