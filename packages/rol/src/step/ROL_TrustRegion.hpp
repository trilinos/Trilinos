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

#include <Teuchos_ScalarTraits.hpp>

namespace ROL { 

enum TrustRegionStepType {
  TrustRegionStep_Secant = 0,
  TrustRegionStep_Newton,
  TrustRegionStep_NewtonKrylov,
  TrustRegionStep_NewtonKrylovSecantPreconditioning
};

enum TrustRegionType {
  TrustRegionType_TruncatedCG = 0,    // TruncatedCG
  TrustRegionType_Dogleg,             // Single Dogleg
  TrustRegionType_DoubleDogleg        // Double Dogleg
};

template<class Real>
class TrustRegion {
private:

  TrustRegionType      type_;
  TrustRegionStepType  step_;

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
  TrustRegion( TrustRegionType type, TrustRegionStepType step,  
              int maxit = 20, Real tol1 = 1.e-4, Real tol2 = 1.e-2, 
              Real delmin = 1.e-8, Real delmax = 5000.0,
              Real eta0 = 0.05, Real eta1 = 0.05, Real eta2 = 0.9,
              Real gamma0 = 0.0625, Real gamma1 = 0.25, Real gamma2 = 2.5, Real TRsafe = 1.0 ) : 
    type_(type), step_(step), maxit_(maxit), tol1_(tol1), tol2_(tol2), 
    delmin_(delmin), delmax_(delmax), eta0_(eta0), eta1_(eta1), eta2_(eta2),
    gamma0_(gamma0), gamma1_(gamma1), gamma2_(gamma2), TRsafe_(TRsafe) {
    eps_ = TRsafe_*Teuchos::ScalarTraits<Real>::eps();
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
      else if (aRed < 0 && pRed_ > 0) {
        flagTR = 2;
      }
      else if (aRed < 0 && pRed_ < 0) { 
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
        if ( secant != Teuchos::null ) {
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
    if ( type_ == TrustRegionType_TruncatedCG ) {
      truncatedCG(s,snorm,del,iflag,iter,x,grad,gnorm,obj,secant);
    }
    //else if ( type_ == TrustRegionType_Dogleg ) {
    //  dogleg(s,snorm,iflag,iter,x,iterTR,grad,gnorm,obj);
    //}
    //else if ( type_ == TrustRegionType_DoubleDogleg ) {
    //  doubledogleg(s,snorm,iflag,iter,x,iterTR,grad,gnorm,obj);
    //}
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
    if ( secant != Teuchos::null && step_ == TrustRegionStep_NewtonKrylovSecantPreconditioning ) { 
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
    if ( secant != Teuchos::null && step_ == TrustRegionStep_Secant ) { 
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

      if ( secant != Teuchos::null && step_ == TrustRegionStep_NewtonKrylovSecantPreconditioning ) { 
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
      if ( secant != Teuchos::null && step_ == TrustRegionStep_Secant ) { 
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

};

}

#endif
