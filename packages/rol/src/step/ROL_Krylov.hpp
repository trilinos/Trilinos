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

#ifndef ROL_KRYLOV_H
#define ROL_KRYLOV_H

/** \class ROL::Krylov
    \brief Provides definitions for Krylov solvers.
*/

#include "ROL_LinearOperator.hpp"
#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

namespace ROL {

template<class Real>
class Krylov {

  bool isInitialized_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > v_;
  Teuchos::RCP<Vector<Real> > p_;
  Teuchos::RCP<Vector<Real> > Hp_;
  Teuchos::RCP<Vector<Real> > MHp_;

  EKrylov ekv_;

  Real tol1_;
  Real tol2_;
  int  maxit_;
  bool useInexact_;

public:
  Krylov( EKrylov ekv = KRYLOV_CG, Real tol1 = 1.e-4, Real tol2 = 1.e-2, 
          int maxit = 100, bool useInexact = false ) 
    : isInitialized_(false), ekv_(ekv), tol1_(tol1), tol2_(tol2), maxit_(maxit), useInexact_(useInexact) {}

  // Run Krylov Method
  void run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
            int &iter, int &flag ) {
    switch ( this->ekv_ ) {
      case KRYLOV_CG: this->CG(x,A,b,M,iter,flag); break;
      case KRYLOV_CR: this->CR(x,A,b,M,iter,flag); break;
      case KRYLOV_LAST: break; // DO NOTHING
    }
  }

  // Use (inexact) CG to solve Newton system 
  void CG( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
           int &iter, int &flag ) {
    if ( !isInitialized_ ) {
      g_  = b.clone();
      v_  = x.clone();
      p_  = x.clone(); 
      Hp_ = b.clone();
      isInitialized_ = true;
    }

    Real rnorm = b.norm(); 
    Real rtol = std::min(tol1_,tol2_*rnorm);
    Real itol = 0.0;

    x.zero(); 
    g_->set(b); 

    itol = 0.0;
    M.apply(*v_, *g_, itol);
    p_->set(*v_); 

    iter = 0; 
    flag = 0;

    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v_->dot(*g_); 

    for (iter = 0; iter < this->maxit_; iter++) {
      if ( this->useInexact_ ) {
        itol = rtol/(this->maxit_ * rnorm); 
      }
      A.apply(*Hp_, *p_, itol);

      kappa = p_->dot(*Hp_);
      if ( kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = gv/kappa;

      x.axpy(alpha,*p_);

      g_->axpy(-alpha,*Hp_);
      rnorm = g_->norm();
      if ( rnorm < rtol ) {
        break;
      }

      itol = 0.0;
      M.apply(*v_, *g_, itol);
      tmp  = gv;         
      gv   = v_->dot(*g_); 
      beta = gv/tmp;
 
      p_->scale(beta);
      p_->axpy(1.0,*v_);
    }
    if ( iter == this->maxit_ ) {
      flag = 1;
    }    
    else {
      iter++;
    }
  }

  // Use (inexact) CR to solve Newton system 
  void CR( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
           int &iter, int &flag ) {
    if ( !isInitialized_ ) {
      g_ = b.clone();
      v_ = x.clone();
      p_ = x.clone();
      Hp_ = b.clone();
      MHp_ = b.clone();
      isInitialized_ = true;
    }

    // Initialize
    Real rnorm = b.norm(); 
    Real rtol = std::min(tol1_,tol2_*rnorm);
    Real itol = 0.0;
    x.zero(); 

    // Apply preconditioner to residual
    itol = 0.0;
    M.apply(*g_,b,itol);

    // Initialize direction p
    p_->set(*g_); 

    // Get Hessian tolerance
    if ( this->useInexact_ ) {
      itol = rtol/(maxit_ * rnorm); 
    }

    // Apply Hessian to residual
    A.apply(*v_, *g_, itol);

    // Apply Hessian to direction p
    //A.apply(*Hp_, *p_, itol);
    Hp_->set(*v_);

    // Initialize scalar quantities
    iter = 0; 
    flag = 0;
    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gHg   = g_->dot(*v_); 

    for (iter = 0; iter < this->maxit_; iter++) {
      itol = 0.0;
      M.apply(*MHp_, *Hp_, itol);
      kappa = Hp_->dot(*MHp_);
      //if ( gHg <= 0.0 || kappa <= 0.0 ) { 
        //flag = 2;
        //break;
      //}
      alpha = gHg/kappa;

      x.axpy(alpha,*p_);

      g_->axpy(-alpha,*MHp_);
      rnorm = g_->norm();
      if ( rnorm < rtol ) {
        break;
      }

      if ( this->useInexact_ ) {
        itol = rtol/(this->maxit_ * rnorm); 
      }
      A.apply(*v_, *g_, itol);
      tmp  = gHg;
      gHg  = g_->dot(*v_);
      beta = gHg/tmp;

      p_->scale(beta);
      p_->axpy(1.0,*g_);

      Hp_->scale(beta);
      Hp_->axpy(1.0,*v_); 
    }
    if ( iter == this->maxit_ ) {
      flag = 1;
    }   
    else {
      iter++;
    } 
  }

  // Run Krylov Method
  void run( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
            ProjectedObjective<Real> &pObj ) {
    switch ( this->ekv_ ) {
      case KRYLOV_CG: this->CG(s,iter,flag,g,x,pObj); break;
      case KRYLOV_CR: this->CR(s,iter,flag,g,x,pObj); break;
      case KRYLOV_LAST: break; // DO NOTHING
    }
  }

  // Use (inexact) CG to solve Newton system 
  void CG( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
           ProjectedObjective<Real> &pObj ) {
    if ( !isInitialized_ ) {
      g_  = g.clone();
      v_  = x.clone();
      p_  = x.clone(); 
      Hp_ = g.clone();
      isInitialized_ = true;
    }

    Real gnorm = g.norm(); 
    Real gtol = std::min(tol1_,tol2_*gnorm);
    Real itol = 0.0;

    s.zero(); 
    g_->set(g); 

    itol = 0.0;
    pObj.reducedPrecond(*v_, *g_, x, g, x, itol);
    p_->set(*v_); 

    iter = 0; 
    flag = 0;

    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v_->dot(*g_); 

    for (iter = 0; iter < this->maxit_; iter++) {
      if ( this->useInexact_ ) {
        itol = gtol/(this->maxit_ * gnorm); 
      }
      pObj.reducedHessVec(*Hp_, *p_, x, g, x, itol);

      kappa = p_->dot(*Hp_);
      if ( kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = gv/kappa;

      s.axpy(alpha,*p_);

      g_->axpy(-alpha,*Hp_);
      gnorm = g_->norm();
      if ( gnorm < gtol ) {
        break;
      }

      itol = 0.0;
      pObj.reducedPrecond(*v_, *g_, x, g, x, itol);
      tmp  = gv;         
      gv   = v_->dot(*g_); 
      beta = gv/tmp;
 
      p_->scale(beta);
      p_->axpy(1.0,*v_);
    }
    if ( iter == this->maxit_ ) {
      flag = 1;
    }    
    else {
      iter++;
    }
  }

  // Use (inexact) CR to solve Newton system 
  void CR( Vector<Real> &s, int &iter, int &flag, const Vector<Real> &g, const Vector<Real> &x, 
           ProjectedObjective<Real> &pObj ) {
    if ( !isInitialized_ ) {
      g_ = g.clone();
      v_ = x.clone();
      p_ = x.clone();
      Hp_ = g.clone();
      MHp_ = g.clone();
      isInitialized_ = true;
    }

    // Initialize
    Real gnorm = g.norm(); 
    Real gtol = std::min(tol1_,tol2_*gnorm);
    Real itol = 0.0;
    s.zero(); 

    // Apply preconditioner to residual
    itol = 0.0;
    pObj.reducedPrecond(*v_, g, x, g, x, itol);

    // Initialize direction p
    p_->set(*v_); 

    // Get Hessian tolerance
    if ( this->useInexact_ ) {
      itol = gtol/(maxit_ * gnorm); 
    }

    // Apply Hessian to residual
    pObj.reducedHessVec(*g_, *v_, x, g, x, itol);

    // Apply Hessian to direction p
    pObj.reducedHessVec(*Hp_, *p_, x, g, x, itol);

    // Initialize scalar quantities
    iter = 0; 
    flag = 0;
    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real vHv   = g_->dot(*v_); 

    for (iter = 0; iter < this->maxit_; iter++) {
      itol = 0.0;
      pObj.reducedPrecond(*MHp_, *Hp_, x, g, x, itol);
      kappa = Hp_->dot(*MHp_);
      //if ( vHv <= 0.0 || kappa <= 0.0 ) { 
      if ( kappa == 0.0 ) {
        flag = 2;
        break;
      }
      alpha = vHv/kappa;

      s.axpy(alpha,*p_);

      v_->axpy(-alpha,*MHp_);
      gnorm = v_->norm();
      if ( gnorm < gtol ) {
        break;
      }

      if ( this->useInexact_ ) {
        itol = gtol/(this->maxit_ * gnorm); 
      }
      pObj.reducedHessVec(*g_, *v_, x, g, x, itol);
      tmp  = vHv;
      vHv  = g_->dot(*v_);
      beta = vHv/tmp;

      p_->scale(beta);
      p_->axpy(1.0,*v_);

      Hp_->scale(beta);
      Hp_->axpy(1.0,*g_); 
    }
    if ( iter == this->maxit_ ) {
      flag = 1;
    }   
    else {
      iter++;
    } 
  }
};


}

#endif
