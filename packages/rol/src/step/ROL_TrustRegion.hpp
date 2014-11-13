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

#ifndef ROL_TRUSTREGION_H
#define ROL_TRUSTREGION_H

/** \class ROL::TrustRegion
    \brief Provides interface for and implements trust-region subproblem solvers.
*/

#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

namespace ROL { 

template<class Real>
class TrustRegion {
private:

  Teuchos::RCP<Vector<Real> > xupdate_;
  Teuchos::RCP<Vector<Real> > s_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > v_;
  Teuchos::RCP<Vector<Real> > p_;
  Teuchos::RCP<Vector<Real> > Hp_;

  ETrustRegion etr_;

  bool useSecantPrecond_;
  bool useSecantHessVec_;

  int maxit_;
  Real tol1_;
  Real tol2_;
  Real delmin_;
  Real delmax_;
  Real eta0_;
  Real eta1_;
  Real eta2_;
  Real gamma0_;
  Real gamma1_;
  Real gamma2_; 

  Real pRed_;

  Real TRsafe_;
  Real eps_;

  Real alpha_;

  std::vector<bool> useInexact_;

  Real ftol_old_;

  Real scale_;
  Real omega_;
  Real force_;
  int  updateIter_;
  int  forceFactor_;
  int cnt_;

  bool useCGTCP_; 

public:

  virtual ~TrustRegion() {}

  // Constructor
  TrustRegion( Teuchos::ParameterList & parlist ) : ftol_old_(ROL_OVERFLOW), cnt_(0), useCGTCP_(false) {
    // Unravel Parameter List
    // Enumerations
    etr_ = StringToETrustRegion( parlist.get("Trust-Region Subproblem Solver Type",  "Cauchy Point"));
    // Trust-Region Parameters
    delmin_ = parlist.get("Minimum Trust-Region Radius",          1.e-8);
    delmax_ = parlist.get("Maximum Trust-Region Radius",          5000.0);
    eta0_   = parlist.get("Step Acceptance Parameter",            0.05);
    eta1_   = parlist.get("Radius Shrinking Threshold",           0.05);
    eta2_   = parlist.get("Radius Growing Threshold",             0.9);
    gamma0_ = parlist.get("Radius Shrinking Rate (Negative rho)", 0.0625);
    gamma1_ = parlist.get("Radius Shrinking Rate (Positive rho)", 0.25);
    gamma2_ = parlist.get("Radius Growing Rate",                  2.5);
    TRsafe_ = parlist.get("Trust-Region Safeguard",               100.0);
    // CG Parameters
    if ( etr_ == TRUSTREGION_TRUNCATEDCG ) {
      maxit_ = parlist.get("Maximum Number of Krylov Iterations",     20);
      tol1_  = parlist.get("Absolute Krylov Tolerance",               1.e-4);
      tol2_  = parlist.get("Relative Krylov Tolerance",               1.e-2);
    }
    eps_    = TRsafe_*ROL_EPSILON;
    alpha_  = -1.0;

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
    scale_       = parlist.get("Value Update Tolerance Scaling",1.e-1);
    omega_       = parlist.get("Value Update Exponent",0.9);
    force_       = parlist.get("Value Update Forcing Sequence Initial Value",1.0);
    updateIter_  = parlist.get("Value Update Forcing Sequence Update Frequency",10);
    forceFactor_ = parlist.get("Value Update Forcing Sequence Reduction Factor",0.1);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &g) {
    xupdate_ = x.clone();
    if (etr_ == TRUSTREGION_TRUNCATEDCG ) {
      s_ = x.clone();
      g_ = g.clone();
      v_ = x.clone();
      p_ = x.clone();
      Hp_ = g.clone();
    }
    else if ( etr_ == TRUSTREGION_DOGLEG ) {
      s_ = x.clone();
      Hp_ = g.clone();
    }
    else if ( etr_ == TRUSTREGION_DOUBLEDOGLEG ) {
      s_ = x.clone();
      v_ = x.clone();
      Hp_ = g.clone();
    }
    else if ( etr_ == TRUSTREGION_CAUCHYPOINT ) {
      Hp_ = g.clone();
      if ( useCGTCP_ ) {
        g_ = g.clone();
        p_ = x.clone();
      }
    }
  }

  void update( Vector<Real> &x, Real &fnew, Real &del, 
               int &nfval, int &ngrad, int &flagTR,
               const Vector<Real> &s, const Real snorm, 
               const Real fold, const Vector<Real> &g, 
               int iter, ProjectedObjective<Real> &pObj ) { 
    Real tol = std::sqrt(ROL_EPSILON);

    // Compute New Function Value
    xupdate_->set(x);
    xupdate_->axpy(1.0,s);
    /***************************************************************************************************/
    // BEGIN INEXACT OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/
    Real fold1 = fold;
    if ( useInexact_[0] ) {
      if ( !(cnt_%updateIter_) && (cnt_ != 0) ) {
        force_ *= forceFactor_;
      }
      Real c = scale_*std::max(1.e-2,std::min(1.0,1.e4*std::max(pRed_,std::sqrt(ROL_EPSILON))));
      Real ftol = c*std::pow(std::min(eta1_,1.0-eta2_)
                   *std::min(std::max(pRed_,std::sqrt(ROL_EPSILON)),force_),1.0/omega_);
      if ( ftol_old_ > ftol || cnt_ == 0 ) {
        ftol_old_ = ftol;
        fold1 = pObj.value(x,ftol_old_);
      }
      pObj.update(*xupdate_,true,iter);
      fnew = pObj.value(*xupdate_,ftol);
      cnt_++;
    }
    else {
      pObj.update(*xupdate_,true,iter);
      fnew = pObj.value(*xupdate_,tol);
    }
    nfval = 1;   
    Real aRed = fold1 - fnew;
    /***************************************************************************************************/
    // FINISH INEXACT OBJECTIVE FUNCTION COMPUTATION
    /***************************************************************************************************/

    // Compute Ratio of Actual and Predicted Reduction
    aRed -= eps_*((1.0 < std::abs(fold1)) ? 1.0 : std::abs(fold1));
    pRed_ -= eps_*((1.0 < std::abs(fold1)) ? 1.0 : std::abs(fold1));
    Real rho  = 0.0; 
    if ((std::abs(aRed) < eps_) && (std::abs(pRed_) < eps_)) {
      rho = 1.0; 
      flagTR = 0;
    }
    else {
      rho = aRed/pRed_;
      if (pRed_ < 0 && aRed > 0) { 
        flagTR = 1;
      }
      else if (aRed <= 0 && pRed_ > 0) {
        flagTR = 2;
      }
      else if (aRed <= 0 && pRed_ < 0) { 
        flagTR = 3;
      }
      else {
        flagTR = 0;
      }
    }

    // Check Sufficient Decrease in the Reduced Quadratic Model
    bool decr = true;
    if ( pObj.isConActivated() && (std::abs(aRed) > eps_) ) { 
      // Compute Criticality Measure || x - P( x - g ) ||
      xupdate_->set(x);
      xupdate_->axpy(-1.0,g.dual());
      pObj.project(*xupdate_);
      xupdate_->scale(-1.0);
      xupdate_->plus(x);
      Real pgnorm = xupdate_->norm();
      // Compute Scaled Measure || x - P( x - lam * PI(g) ) ||
      xupdate_->set(g.dual());
      pObj.pruneActive(*xupdate_,g,x);
      Real lam = std::min(1.0, del/xupdate_->norm());
      xupdate_->scale(-lam);
      xupdate_->plus(x);
      pObj.project(*xupdate_);
      xupdate_->scale(-1.0);
      xupdate_->plus(x);      
      pgnorm *= xupdate_->norm();
      // Sufficient decrease?
      decr = ( aRed >= 0.1*eta0_*pgnorm );
    }
    
    // Accept or Reject Step and Update Trust Region
    if ((rho < eta0_ && flagTR == 0) || flagTR >= 2 || !decr ) { // Step Rejected 
      pObj.update(x,true,iter);
      fnew = fold1;
      if (rho < 0.0) { // Negative reduction, interpolate to find new trust-region radius
        Real gs = s.dot(g.dual());
        pObj.hessVec(*Hp_,s,x,tol);
        Real modelVal = s.dot(Hp_->dual());
        modelVal *= 0.5;
        modelVal += gs + fold1;
        Real theta = (1.0-eta2_)*gs/((1.0-eta2_)*(fold1+gs)+eta2_*modelVal-fnew);
        del = std::min(gamma1_*snorm,std::max(gamma0_,theta)*del);
      }
      else { // Shrink trust-region radius
        del = gamma1_*snorm; 
      } 
    }
    else if ((rho >= eta0_ && flagTR != 3) || flagTR == 1) { // Step Accepted
      x.axpy(1.0,s);
      if (rho >= eta2_) { // Increase trust-region radius
        del = std::min(gamma2_*del,delmax_);
      }
    }
  }

  void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
            const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) { 
    // Run Trust Region
    if ( etr_ == TRUSTREGION_CAUCHYPOINT ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
    }
    else if ( etr_ == TRUSTREGION_TRUNCATEDCG ) {
      truncatedCG(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
      if ( pObj.isConActivated() ) {
        pObj.computeProjectedStep(s,x);
      }
#if 0
      if ( pObj.isConActivated() ) {
        truncatedCG_proj(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
      }
      else {
        truncatedCG(s,snorm,del,iflag,iter,x,grad,grad,gnorm,pObj);
      }
#endif
    }
    else if ( etr_ == TRUSTREGION_DOGLEG ) {
      dogleg(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
    }
    else if ( etr_ == TRUSTREGION_DOUBLEDOGLEG ) {
      doubledogleg(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
    }
#if 0
    // If constraints are active, then determine a feasible step
    if ( pObj.isConActivated() && etr_ != TRUSTREGION_CAUCHYPOINT ) {
      // Compute projected step stmp = P( x + alpha*s ) - x and xtmp = x + stmp
      Real alpha = 1.0;
      Teuchos::RCP<Vector<Real> > stmp = x.clone();
      stmp->set(s);
      stmp->scale(alpha);
      pObj.computeProjectedStep(*stmp,x);
      Teuchos::RCP<Vector<Real> > xtmp = x.clone();
      xtmp->set(x);
      xtmp->axpy(1.0,*stmp);
      // Compute model components for alpha = 1.0
      Real tol   = std::sqrt(ROL_EPSILON);
      Teuchos::RCP<Vector<Real> > Bs = x.clone();
      pObj.hessVec(*Bs,*stmp,x,tol);
      Real sBs   = Bs->dot(*stmp);
      Real gs    = grad.dot(*stmp);
      Real val   = gs + 0.5*sBs;
      Real val0  = val;
      // Backtrack alpha until x+alpha*s is feasible
      int cnt   = 0;
      int maxit = 10;
      while ( val > val0 || !pObj.isFeasible(*xtmp) ) { 
        // Backtrack alpha
        alpha *= 0.5;
        // Computed projected step P( x + alpha*s ) - x and xtmp = x + stmp
        stmp->set(s);
        stmp->scale(alpha);
        pObj.computeProjectedStep(*stmp,x);        
        xtmp->set(x);
        xtmp->axpy(1.0,*stmp);
        // Evaluate Model
        val0 = val;
        pObj.hessVec(*Bs,*stmp,x,tol);
        sBs = Bs->dot(*stmp);
        gs  = grad.dot(*stmp);
        val = gs + 0.5*sBs;
        // Update iteration count
        cnt++;
        if ( cnt >= maxit ) { break; }
      }
      s.set(*stmp);
      pRed_ = -val;
    }
#endif
  }

  void cauchypoint( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
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
  }

  void truncatedCG( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
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

  void dogleg( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
               const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Compute quasi-Newton step
    pObj.reducedInvHessVec(*s_,grad,x,grad,x,tol);
    s_->scale(-1.0);
    Real sNnorm = s_->norm();
    Real gsN    = s_->dot(grad.dual());
    bool negCurv = false;
    if ( gsN >= 0.0 ) {
      negCurv = true;
    }

    if ( negCurv ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
      iflag = 2;
    }  
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) {        // Use the quasi-Newton step
        s.set(*s_); 
        snorm = sNnorm;
        pRed_ = -0.5*gsN;
        iflag = 0;
      }
      else {                      // quasi-Newton step is outside of trust region
        pObj.reducedHessVec(*Hp_,grad,x,grad,x,tol);
        Real alpha  = 0.0;
        Real beta   = 0.0;
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = grad.dot(*Hp_);
        Real gamma  = gnorm2/gBg;
        if ( gamma*gnorm >= del || gBg <= 0.0 ) {
            alpha = 0.0;
            beta  = del/gnorm;
            s.set(grad.dual()); 
            s.scale(-beta); 
            snorm = del;
            iflag = 2;
        }
        else {
          Real a = sNnorm*sNnorm + 2.0*gamma*gsN + gamma*gamma*gnorm2;
          Real b = -gamma*gsN - gamma*gamma*gnorm2;
          Real c = gamma*gamma*gnorm2 - del*del;
          alpha  = (-b + sqrt(b*b - a*c))/a;
          beta   = gamma*(1.0-alpha);
          s.set(grad.dual());
          s.scale(-beta);
          s.axpy(alpha,*s_);
          snorm = del;
          iflag = 1;
        }
        pRed_ = (alpha*(0.5*alpha-1)*gsN - 0.5*beta*beta*gBg + beta*(1-alpha)*gnorm2);
      }
    }
  }

  void doubledogleg( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                     const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Compute quasi-Newton step
    pObj.reducedInvHessVec(*s_,grad,x,grad,x,tol);
    s_->scale(-1.0);
    Real sNnorm = s_->norm();
    Real tmp    = s_->dot(grad.dual());
    bool negCurv = false;
    if ( tmp >= 0.0 ) {
      negCurv = true;
    }
    Real gsN = std::abs(tmp);

    if ( negCurv ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,pObj);
      iflag = 2;
    }  
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) {        // Use the quasi-Newton step
        s.set(*s_); 
        snorm = sNnorm;
        pRed_ = 0.5*gsN;
        iflag = 0;
      }
      else {                      // quasi-Newton step is outside of trust region
        pObj.reducedHessVec(*Hp_,grad,x,grad,x,tol);
        Real alpha  = 0.0;
        Real beta   = 0.0;
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = grad.dot(*Hp_);
        Real gamma1 = gnorm/gBg;
        Real gamma2 = gnorm/gsN;
        Real eta    = 0.8*gamma1*gamma2 + 0.2;
        if (eta*sNnorm <= del || gBg <= 0.0) { // Dogleg Point is inside trust region
          alpha = del/sNnorm;
          beta  = 0.0;
          s.set(*s_);
          s.scale(alpha);
          snorm = del;
          iflag = 1;
        }
        else {
          if (gnorm2*gamma1 >= del) { // Cauchy Point is outside trust region
            alpha = 0.0;
            beta  = -del/gnorm;
            s.set(grad.dual()); 
            s.scale(beta); 
            snorm = del;
            iflag = 2;
          }
          else {              // Find convex combination of Cauchy and Dogleg point
            s.set(grad.dual());
            s.scale(-gamma1*gnorm);
            v_->set(s);
            v_->scale(-1.0);
            v_->axpy(eta,*s_);
            Real wNorm = v_->dot(*v_);
            Real sigma = del*del-std::pow(gamma1*gnorm,2.0);
            Real phi   = s.dot(*v_);
            Real theta = (-phi + std::sqrt(phi*phi+wNorm*sigma))/wNorm;
            s.axpy(theta,*v_); 
            snorm = del;
            alpha = theta*eta;
            beta  = (1.0-theta)*(-gamma1*gnorm);
            iflag = 3;
          }
        }
        pRed_ = -(alpha*(0.5*alpha-1)*gsN + 0.5*beta*beta*gBg + beta*(1-alpha)*gnorm2);
      }
    }
  }

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
