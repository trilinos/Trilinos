// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUNCATEDCG_H
#define ROL_TRUNCATEDCG_H

/** \class ROL::TruncatedCG
    \brief Provides interface for truncated CG trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class TruncatedCG : public TrustRegion<Real> {
private:
  ROL::Ptr<Vector<Real> > primalVector_;

  ROL::Ptr<Vector<Real> > s_;
  ROL::Ptr<Vector<Real> > g_;
  ROL::Ptr<Vector<Real> > v_;
  ROL::Ptr<Vector<Real> > p_;
  ROL::Ptr<Vector<Real> > Hp_;

  int maxit_;
  Real tol1_;
  Real tol2_;

  Real pRed_;

public:

  // Constructor
  TruncatedCG( ROL::ParameterList &parlist ) : TrustRegion<Real>(parlist), pRed_(0) {
    // Unravel Parameter List
    Real em4(1e-4), em2(1e-2);
    maxit_ = parlist.sublist("General").sublist("Krylov").get("Iteration Limit",20);
    tol1_  = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance",em4);
    tol2_  = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance",em2);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);

    primalVector_ = x.clone();

    s_ = s.clone();
    g_ = g.clone();
    v_ = s.clone();
    p_ = s.clone();
    Hp_ = g.clone();
  }

  void run( Vector<Real>           &s,
            Real                   &snorm,
            int                    &iflag,
            int                    &iter,
            const Real              del,
            TrustRegionModel<Real> &model ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    const Real zero(0), one(1), two(2), half(0.5);
    // Initialize step
    s.zero(); s_->zero();
    snorm = zero;
    Real snorm2(0), s1norm2(0);
    // Compute (projected) gradient
    model.dualTransform(*g_,*model.getGradient());
    Real gnorm = g_->norm(), normg = gnorm;
    const Real gtol = std::min(tol1_,tol2_*gnorm);
    // Preconditioned (projected) gradient vector
    model.precond(*v_,*g_,s,tol);
    // Initialize basis vector
    p_->set(*v_); p_->scale(-one);
    Real pnorm2 = v_->dot(g_->dual());
    if ( pnorm2 <= zero ) {
      iflag = 4;
      iter  = 0;
      return;
    }
    // Initialize scalar storage
    iter = 0; iflag = 0;
    Real kappa(0), beta(0), sigma(0), alpha(0), tmp(0), sMp(0);
    Real gv = v_->dot(g_->dual());
    pRed_ = zero;
    // Iterate CG
    for (iter = 0; iter < maxit_; iter++) {
      // Apply Hessian to direction p
      model.hessVec(*Hp_,*p_,s,tol);
      // Check positivity of Hessian
      kappa = p_->dot(Hp_->dual());
      if (kappa <= zero) {
        sigma = (-sMp+sqrt(sMp*sMp+pnorm2*(del*del-snorm2)))/pnorm2;
        s.axpy(sigma,*p_);
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
        iflag = 3; 
        break;
      }
      // Update model predicted reduction
      pRed_ += half*alpha*gv;
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
      gv    = v_->dot(g_->dual());
      beta  = gv/tmp;    
      // Update basis vector
      p_->scale(beta);
      p_->axpy(-one,*v_);
      sMp    = beta*(sMp+alpha*pnorm2);
      pnorm2 = gv + beta*beta*pnorm2; 
    }
    // Update model predicted reduction
    if (iflag > 0) {
      pRed_ += sigma*(gv-half*sigma*kappa);
    }
    // Check iteration count
    if (iter == maxit_) {
      iflag = 1;
    }
    if (iflag != 1) { 
      iter++;
    }
    // Update norm of step and update model predicted reduction
    model.primalTransform(*s_,s);
    s.set(*s_);
    snorm = s.norm();
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }

#if 0
  void truncatedCG_proj( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                         const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());

    const Real gtol = std::min(tol1_,tol2_*gnorm);

    // Compute Cauchy Point
    Real scnorm = 0.0;
    ROL::Ptr<Vector<Real> > sc = x.clone();
    cauchypoint(*sc,scnorm,del,iflag,iter,x,grad,gnorm,pObj);
    ROL::Ptr<Vector<Real> > xc = x.clone();
    xc->set(x);
    xc->plus(*sc);

    // Old and New Step Vectors
    s.set(*sc); 
    snorm = s.norm();
    Real snorm2  = snorm*snorm;
    ROL::Ptr<Vector<Real> > s1 = x.clone();
    s1->zero();
    Real s1norm2 = 0.0;

    // Gradient Vector
    ROL::Ptr<Vector<Real> > g = x.clone(); 
    g->set(grad);
    ROL::Ptr<Vector<Real> > Hs = x.clone();
    pObj.reducedHessVec(*Hs,s,*xc,x,tol);
    g->plus(*Hs);
    Real normg = g->norm();

    // Preconditioned Gradient Vector
    ROL::Ptr<Vector<Real> > v  = x.clone();
    pObj.reducedPrecond(*v,*g,*xc,x,tol);

    // Basis Vector
    ROL::Ptr<Vector<Real> > p = x.clone(); 
    p->set(*v); 
    p->scale(-1.0);
    Real pnorm2 = v->dot(*g);

    // Hessian Times Basis Vector
    ROL::Ptr<Vector<Real> > Hp = x.clone();

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
