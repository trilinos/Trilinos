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

#include <Teuchos_ScalarTraits.hpp>

namespace ROL {

template<class Real>
class Krylov {

  bool useSecant_;
  Teuchos::RCP<Secant<Real> > secant_;

  Real eps_;
  Real tol1_;
  Real tol2_;
  int  maxit_;

public:
  Krylov( Real tol1 = 1.e-4, Real tol2 = 1.e-2, int maxit = 100 ) : tol1_(tol1), tol2_(tol2), maxit_(maxit) {
    eps_ = Teuchos::ScalarTraits<Real>::eps();
    useSecant_ = false;
  }
  Krylov( SecantType type, Real tol1 = 1.e-4, Real tol2 = 1.e-2, int maxit = 100, int L = 10, int BBtype = 1 ) 
    : tol1_(tol1), tol2_(tol2), maxit_(maxit) {
    eps_ = Teuchos::ScalarTraits<Real>::eps();

    useSecant_ = true;
    if ( type == Secant_lBFGS ) {
      secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
    }
    else if ( type == Secant_lDFP ) {
      secant_ = Teuchos::rcp( new lDFP<Real>(L) );
    }
    else if ( type == Secant_BarzilaiBorwein ) {
      secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
    }
  }

  // Use CG to solve Newton system
  void CG( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, Objective<Real> &obj, int it = 0 ) {
    Real gtol = std::min(tol1_,tol2_*g.norm());

    s.zero(); 

    Teuchos::RCP<Vector<Real> > gnew = x.clone(); 
    gnew->set(g); 

    Teuchos::RCP<Vector<Real> > v = x.clone();  
    if ( useSecant_ ) {
      secant_->applyB( *v, *gnew, x, it );
    }
    else {
      obj.precond( *v, *gnew, x );  
    }

    Teuchos::RCP<Vector<Real> > p = x.clone(); 
    p->set(*v); 

    Teuchos::RCP<Vector<Real> > Hp = x.clone();
    obj.hessVec( *Hp, *p, x );  

    iter = 0; 
    flag = 0;

    Real gnorm = 0;
    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v->dot(*gnew); 

    for (iter = 0; iter < maxit_; iter++) {
      kappa = p->dot(*Hp);
      if ( kappa <= eps_ ) { 
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

      if ( useSecant_ ) {
        secant_->applyB( *v, *gnew, x, it );
      }
      else {
        obj.precond( *v, *gnew, x );  
      }
      tmp  = gv;         
      gv   = v->dot(*gnew); 
      beta = gv/tmp;  
 
      p->scale(beta);
      p->axpy(1.0,*v);

      obj.hessVec( *Hp, *p, x );
    }
    iter++;
    if ( iter == maxit_ ) {
      flag = 1;
    }    
  }

};

}

#endif
