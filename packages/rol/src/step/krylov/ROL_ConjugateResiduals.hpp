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

#ifndef ROL_CONJUGATERESIDUALS_H
#define ROL_CONJUGATERESIDUALS_H

/** \class ROL::ConjugateResiduals
    \brief Provides definition of the Conjugate Residual solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class ConjugateResiduals : public Krylov<Real> {

  bool isInitialized_;
  bool useInexact_;
  ROL::Ptr<Vector<Real> > r_;
  ROL::Ptr<Vector<Real> > v_;
  ROL::Ptr<Vector<Real> > p_;
  ROL::Ptr<Vector<Real> > Ap_;
  ROL::Ptr<Vector<Real> > MAp_;

public:
  ConjugateResiduals( Real absTol = 1.e-4, Real relTol = 1.e-2, int maxit = 100, bool useInexact = false ) 
    : Krylov<Real>(absTol,relTol,maxit), isInitialized_(false), useInexact_(useInexact) {}

  // Run Krylov Method
  Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
            int &iter, int &flag ) {
    if ( !isInitialized_ ) {
      r_   = x.clone();
      v_   = b.clone();
      p_   = x.clone();
      Ap_  = b.clone();
      MAp_ = x.clone();
      isInitialized_ = true;
    }

    // Initialize
    Real rnorm = b.norm(); 
    Real rtol = std::min(Krylov<Real>::getAbsoluteTolerance(),Krylov<Real>::getRelativeTolerance()*rnorm);
    Real itol = std::sqrt(ROL_EPSILON<Real>());
    x.zero(); 

    // Apply preconditioner to residual
    M.applyInverse(*r_,b,itol);

    // Initialize direction p
    p_->set(*r_); 

    // Get Hessian tolerance
    if ( useInexact_ ) {
      itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm); 
    }

    // Apply Hessian to residual
    A.apply(*v_, *r_, itol);

    // Apply Hessian to direction p
    //A.apply(*Ap_, *p_, itol);
    Ap_->set(*v_);

    // Initialize scalar quantities
    iter = 0; 
    flag = 0;
    Real kappa(0), beta(0), alpha(0), tmp(0);
    Real gHg   = r_->dot(v_->dual()); 

    for (iter = 0; iter < (int)Krylov<Real>::getMaximumIteration(); iter++) {
      itol = std::sqrt(ROL_EPSILON<Real>());
      M.applyInverse(*MAp_, *Ap_, itol);
      kappa = MAp_->dot(Ap_->dual());
      //if ( gHg <= 0.0 || kappa <= 0.0 ) { 
        //flag = 2;
        //break;
      //}
      alpha = gHg/kappa;

      x.axpy(alpha,*p_);

      r_->axpy(-alpha,*MAp_);
      rnorm = r_->norm();
      if ( rnorm < rtol ) {
        break;
      }

      if ( useInexact_ ) {
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm); 
      }
      A.apply(*v_, *r_, itol);
      tmp  = gHg;
      gHg  = r_->dot(v_->dual());
      beta = gHg/tmp;

      p_->scale(beta);
      p_->plus(*r_);

      Ap_->scale(beta);
      Ap_->plus(*v_); 
    }
    if ( iter == (int)Krylov<Real>::getMaximumIteration() ) {
      flag = 1;
    }   
    else {
      iter++;
    } 
    return rnorm;
  }
};


}

#endif
