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

#ifndef ROL_CONJUGATEGRADIENTS_H
#define ROL_CONJUGATEGRADIENTS_H

/** \class ROL::ConjugateGradients
    \brief Provides definitions of the Conjugate Gradient solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class ConjugateGradients : public Krylov<Real> {

  bool isInitialized_;
  bool useInexact_;
  Teuchos::RCP<Vector<Real> > r_;
  Teuchos::RCP<Vector<Real> > v_;
  Teuchos::RCP<Vector<Real> > p_;
  Teuchos::RCP<Vector<Real> > Ap_;

public:
  ConjugateGradients(Real absTol = 1.e-4, Real relTol = 1.e-2, unsigned maxit = 100, bool useInexact = false) 
    : Krylov<Real>(absTol,relTol,maxit), isInitialized_(false), useInexact_(useInexact) {}

  void run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
            int &iter, int &flag ) {
    if ( !isInitialized_ ) {
      r_  = b.clone();
      v_  = x.clone();
      p_  = x.clone(); 
      Ap_ = b.clone();
      isInitialized_ = true;
    }

    Real rnorm = b.norm(); 
    Real rtol = std::min(Krylov<Real>::getAbsoluteTolerance(),Krylov<Real>::getRelativeTolerance()*rnorm);
    Real itol = 0.0;

    x.zero(); 
    r_->set(b); 

    itol = 0.0;
    M.apply(*v_, *r_, itol);
    p_->set(*v_); 

    iter = 0; 
    flag = 0;

    Real kappa = 0.0; 
    Real beta  = 0.0; 
    Real alpha = 0.0; 
    Real tmp   = 0.0;
    Real gv    = v_->dot(r_->dual()); 

    for (iter = 0; iter < (int)Krylov<Real>::getMaximumIteration(); iter++) {
      if ( useInexact_ ) {
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm); 
      }
      A.apply(*Ap_, *p_, itol);

      kappa = p_->dot(Ap_->dual());
      if ( kappa <= 0.0 ) { 
        flag = 2;
        break;
      }
      alpha = gv/kappa;

      x.axpy(alpha,*p_);

      r_->axpy(-alpha,*Ap_);
      rnorm = r_->norm();
      if ( rnorm < rtol ) {
        break;
      }

      itol = 0.0;
      M.apply(*v_, *r_, itol);
      tmp  = gv;         
      gv   = v_->dot(r_->dual()); 
      beta = gv/tmp;
 
      p_->scale(beta);
      p_->axpy(1.0,*v_);
    }
    if ( iter == (int)Krylov<Real>::getMaximumIteration() ) {
      flag = 1;
    }    
    else {
      iter++;
    }
  }
};


}

#endif
