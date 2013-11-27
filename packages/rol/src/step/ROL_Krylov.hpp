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

#ifndef ROL_KRYLOV_H
#define ROL_KRYLOV_H

/** \class ROL::Kryolv
    \brief Provides defintions for Krylov solvers.
*/

namespace ROL {

template<class Real>
class Krylov {

  Real tol1_;
  Real tol2_;
  int  maxit_;
  bool useInexact_;

public:
  Krylov( Real tol1 = 1.e-4, Real tol2 = 1.e-2, int maxit = 100, bool useInexact = false ) 
    : tol1_(tol1), tol2_(tol2), maxit_(maxit), useInexact_(useInexact) {}

  // Use (inexact) CG to solve Newton system 
  void CG( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
           Objective<Real> &obj, Teuchos::RCP<Secant<Real> > secant = Teuchos::null ) {
    Real gnorm = g.norm(); 
    Real gtol = std::min(tol1_,tol2_*gnorm);
    Real itol = 0.0;

    s.zero(); 

    Teuchos::RCP<Vector<Real> > gnew = x.clone(); 
    gnew->set(g); 

    Teuchos::RCP<Vector<Real> > v = x.clone();  
    if ( secant != Teuchos::null ) {
      secant->applyB( *v, *gnew, x );
    }
    else {
      obj.precond( *v, *gnew, x );  
    }

    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 

    Teuchos::RCP<Vector<Real> > Hp = x.clone();
    itol = 0.0;
    if ( useInexact_ ) {
      itol = gtol/(maxit_ * gnorm); 
    }
    obj.hessVec( *Hp, *p, x, itol );  

    iter = 0; 
    flag = 0;

    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v->dot(*gnew); 

    for (iter = 0; iter < maxit_; iter++) {
      kappa = p->dot(*Hp);
      if ( kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = gv/kappa;

      s.axpy(alpha,*p);

      gnew->axpy(-alpha,*Hp);
      gnorm = gnew->norm();
      if ( gnorm < gtol ) {
        break;
      }

      if ( secant != Teuchos::null ) {
        secant->applyB( *v, *gnew, x );
      }
      else {
        obj.precond( *v, *gnew, x );  
      }
      tmp  = gv;         
      gv   = v->dot(*gnew); 
      beta = gv/tmp;
 
      p->scale(beta);
      p->axpy(1.0,*v);

      itol = 0.0;
      if ( useInexact_ ) {
        itol = gtol/(maxit_ * gnorm); 
      }
      obj.hessVec( *Hp, *p, x, itol );
    }
    iter++;
    if ( iter == maxit_ ) {
      flag = 1;
    }    
  }
};

}

#endif
