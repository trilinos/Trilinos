//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_TRUSTREGION_H
#define ROL_TRUSTREGION_H

/** \class ROL::TrustRegion
    \brief Provides interfrace for different trust-region subproblem solvers.
*/

#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class TrustRegion {
private:

  ETrustRegion etr_;
  EDescent     edesc_;

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

public:

  virtual ~TrustRegion() {}

  // Constructor
  TrustRegion( ETrustRegion etr, EDescent edesc,  
              int maxit = 20, Real tol1 = 1.e-4, Real tol2 = 1.e-2, 
              Real delmin = 1.e-8, Real delmax = 5000.0,
              Real eta0 = 0.05, Real eta1 = 0.05, Real eta2 = 0.9,
              Real gamma0 = 0.0625, Real gamma1 = 0.25, Real gamma2 = 2.5, Real TRsafe = 1.0 ) : 
    etr_(etr), edesc_(edesc), maxit_(maxit), tol1_(tol1), tol2_(tol2), 
    delmin_(delmin), delmax_(delmax), eta0_(eta0), eta1_(eta1), eta2_(eta2),
    gamma0_(gamma0), gamma1_(gamma1), gamma2_(gamma2), TRsafe_(TRsafe) {
    eps_ = TRsafe_*ROL_EPSILON;
  }

  void update( Vector<Real> &x, Real &fnew, Real &del, 
               int &nfval, int &ngrad, int &flagTR,
               const Vector<Real> &s, const Real snorm, 
               const Real fold, const Vector<Real> &g, Objective<Real> &obj,
               Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) { 
    // Compute New Function Value
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(1.0,s);
    fnew = obj.value(*xnew);
    nfval = 1;   

    // Compute Ratio of Actual and Predicted Reduction
    Real aRed = fold - fnew;
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

    // Accept or Reject Step and Update Trust Region
    if ((rho < eta0_ && flagTR == 0) || flagTR >= 2 ) {      // Step Rejected 
      fnew = fold;
      if (rho < 0.0) {  
        Real gs = g.dot(s);
        Teuchos::RCP<Vector<Real> > Hs = x.clone();
        if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) {
          secant->applyB(*Hs,s,x);
        } 
        else {
          obj.hessVec(*Hs,s,x);
        }
        Real modelVal = Hs->dot(s);
        modelVal *= 0.5;
        modelVal += gs;
        Real theta = (1.0-eta2_)*gs/((1.0-eta2_)*(fold+gs)+eta2_*modelVal-fnew);
        del  = std::min(gamma1_*snorm,std::max(gamma0_,theta)*del);
      }
      else { 
        del = gamma1_*snorm; 
      } 
    }
    else if ((rho >= eta0_ && flagTR != 3) || flagTR == 1) { // Step Accepted
      x.axpy(1.0,s);
    }
    if ((rho >= eta0_ && flagTR != 3) || flagTR == 1) {
      if (rho >= eta2_) {
        del = std::min(gamma2_*del,delmax_);
      }
    }
  }

  void run( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
            const Vector<Real> &grad, const Real &gnorm, Objective<Real> &obj, 
            Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) { 
    // Run Trust Region
    if ( edesc_ == DESCENT_STEEPEST || etr_ == TRUSTREGION_CAUCHYPOINT ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
    }
    else {
      if ( etr_ == TRUSTREGION_TRUNCATEDCG ) {
        truncatedCG(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
      }
      else if ( etr_ == TRUSTREGION_DOGLEG ) {
        dogleg(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
      }
      else if ( etr_ == TRUSTREGION_DOUBLEDOGLEG ) {
        doubledogleg(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
      }
    }
  }

  void cauchypoint( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                    const Vector<Real> &grad, const Real &gnorm, Objective<Real> &obj,
                    Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) {
    Teuchos::RCP<Vector<Real> > Bg = x.clone();
    if ( secant != Teuchos::null ) {
      secant->applyB(*Bg,grad,x);
    }
    else {
      obj.hessVec(*Bg,grad,x);
    }
    Real gBg = Bg->dot(grad);
    Real tau = 1.0;
    if ( gBg > 0.0 ) {
      tau = std::min(1.0, gnorm*gnorm*gnorm/gBg);
    }

    s.set(grad);
    s.scale(-tau*del/gnorm);
    snorm = tau*del;
    iflag = 0;
    iter  = 0;
    pRed_ = tau*del/gnorm * pow(gnorm,2.0) - 0.5*pow(tau*del/gnorm,2.0)*gBg;
  }

  void truncatedCG( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                    const Vector<Real> &grad, const Real &gnorm, Objective<Real> &obj,
                    Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) {

    const Real gtol = std::min(tol2_,tol1_*gnorm);

    // Old and New Step Vectors
    s.zero(); 
    snorm = 0.0;
    Real snorm2 = 0.0;
    Teuchos::RCP<Vector<Real> > s1 = x.clone();
    s1->zero();
    Real s1norm2 = 0.0;

    // Gradient Vector
    Teuchos::RCP<Vector<Real> > g = x.clone(); 
    g->set(grad); 
    Real normg = gnorm;

    // Preconditioned Gradient Vector
    Teuchos::RCP<Vector<Real> > v  = x.clone();
    if ( secant != Teuchos::null && edesc_ == DESCENT_SECANTPRECOND ) { 
      secant->applyH(*v,*g,x); 
    }
    else { 
      obj.precond( *v, *g, x );        
    }

    // Basis Vector
    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 
    p->scale(-1.0);
    Real pnorm2 = v->dot(*g);

    // Hessian Times Basis Vector
    Teuchos::RCP<Vector<Real> > Hp = x.clone();
    if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) { 
      secant->applyB(*Hp,*p,x); 
    }
    else { 
      obj.hessVec( *Hp, *p, x );
    } 

    iter       = 0; 
    iflag      = 0;
    Real kappa = 0.0;
    Real beta  = 0.0; 
    Real sigma = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v->dot(*g);
    Real sMp   = 0.0;
    pRed_      = 0.0;

    for (iter = 0; iter < maxit_; iter++) {
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
        iflag = 1; 
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

      if ( secant != Teuchos::null && edesc_ == DESCENT_SECANTPRECOND ) { 
        secant->applyH(*v,*g,x); 
      }
      else { 
        obj.precond( *v, *g, x );        
      }
      tmp   = gv; 
      gv    = v->dot(*g);
      beta  = gv/tmp;    

      p->scale(beta);
      p->axpy(-1.0,*v);
      sMp    = beta*(sMp+alpha*pnorm2);
      pnorm2 = gv + beta*beta*pnorm2; 
      if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) { 
        secant->applyB(*Hp,*p,x); 
      }
      else { 
        obj.hessVec( *Hp, *p, x );         
      }
    }
    if (iflag > 0) {
      pRed_ += sigma*(gv-0.5*sigma*kappa);
    }

    if (iter == maxit_) {
      iflag = 3;
    }
    if (iflag != 3) { 
      iter++;
    }

    snorm = s.norm();
  }

  void dogleg( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
               const Vector<Real> &grad, const Real &gnorm, Objective<Real> &obj,
               Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) {
    // Compute quasi-Newton step
    Teuchos::RCP<Vector<Real> > sN = x.clone();
    if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) {
      secant->applyH(*sN,grad,x); 
    }
    else {
      obj.invHessVec(*sN,grad,x);
    }
    sN->scale(-1.0);
    Real sNnorm = sN->norm();
    Real gsN    = grad.dot(*sN);
    bool negCurv = false;
    if ( gsN >= 0.0 ) {
      negCurv = true;
    }

    if ( negCurv ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
      iflag = 2;
    }  
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) {        // Use the quasi-Newton step
        s.set(*sN); 
        snorm = sNnorm;
        pRed_ = -0.5*gsN;
        iflag = 0;
      }
      else {                      // quasi-Newton step is outside of trust region
        Teuchos::RCP<Vector<Real> > Bg = x.clone(); 
        if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) {
          secant->applyB(*Bg,grad,x);
        }
        else {
          obj.hessVec(*Bg,grad,x);
        }
        Real alpha  = 0.0;
        Real beta   = 0.0;
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = grad.dot(*Bg);
        Real gamma  = gnorm2/gBg;
        if ( gamma*gnorm >= del || gBg <= 0.0 ) {
            alpha = 0.0;
            beta  = del/gnorm;
            s.set(grad); 
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
          s.set(grad);
          s.scale(-beta);
          s.axpy(alpha,*sN);
          snorm = del;
          iflag = 1;
        }
        pRed_ = (alpha*(0.5*alpha-1)*gsN - 0.5*beta*beta*gBg + beta*(1-alpha)*gnorm2);
      }
    }
  }

  void doubledogleg( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
                     const Vector<Real> &grad, const Real &gnorm, Objective<Real> &obj,
                     Teuchos::RCP<Secant<Real> > &secant = Teuchos::null ) {
    // Compute quasi-Newton step
    Teuchos::RCP<Vector<Real> > sN = x.clone();
    if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) {
      secant->applyH(*sN,grad,x); 
    }
    else {
      obj.invHessVec(*sN,grad,x);
    }
    sN->scale(-1.0);
    Real sNnorm = sN->norm();
    Real tmp    = grad.dot(*sN);
    bool negCurv = false;
    if ( tmp >= 0.0 ) {
      negCurv = true;
    }
    Real gsN = std::abs(tmp);

    if ( negCurv ) {
      cauchypoint(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
    }  
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) {        // Use the quasi-Newton step
        s.set(*sN); 
        snorm = sNnorm;
        pRed_ = 0.5*gsN;
        iflag = 0;
      }
      else {                      // quasi-Newton step is outside of trust region
        Teuchos::RCP<Vector<Real> > Bg = x.clone(); 
        if ( secant != Teuchos::null && edesc_ == DESCENT_SECANT ) {
          secant->applyB(*Bg,grad,x);
        }
        else {
          obj.hessVec(*Bg,grad,x);
        }
        Real alpha  = 0.0;
        Real beta   = 0.0;
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = grad.dot(*Bg);
        Real gamma1 = gnorm/gBg;
        Real gamma2 = gnorm/gsN;
        Real eta    = 0.8*gamma1*gamma2 + 0.2;
        if (eta*sNnorm <= del || gBg <= 0.0) { // Dogleg Point is inside trust region
          alpha = del/sNnorm;
          beta  = 0.0;
          s.set(*sN);
          s.scale(alpha);
          snorm = del;
          iflag = 1;
        }
        else {
          if (gnorm2*gamma1 >= del) { // Cauchy Point is outside trust region
            alpha = 0.0;
            beta  = -del/gnorm;
            s.set(grad); 
            s.scale(beta); 
            snorm = del;
            iflag = 2;
          }
          else {              // Find convex combination of Cauchy and Dogleg point
            s.set(grad);
            s.scale(-gamma1*gnorm);
            Teuchos::RCP<Vector<Real> > w = x.clone(); 
            w->set(s);
            w->scale(-1.0);
            w->axpy(eta,*sN);
            Real wNorm = w->dot(*w);
            Real sigma = del*del-std::pow(gamma1*gnorm,2.0);
            Real phi   = s.dot(*w);
            Real theta = (-phi + std::sqrt(phi*phi+wNorm*sigma))/wNorm;
            s.axpy(theta,*w); 
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

};

}

#endif
