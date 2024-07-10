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

#ifndef ROL_BICGSTAB_H
#define ROL_BICGSTAB_H

/** \class ROL::ConjugateGradients
    \brief Provides definitions of the Conjugate Gradient solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class BiCGSTAB : public Krylov<Real> {

  bool isInitialized_;
  const bool useInexact_;
  Ptr<Vector<Real>> r_, r1_, p_, v_, s_, t_, h_, y_, z_;

public:
  BiCGSTAB(Real absTol = 1.e-4, Real relTol = 1.e-2, unsigned maxit = 100, bool useInexact = false)
    : Krylov<Real>(absTol,relTol,maxit), isInitialized_(false), useInexact_(useInexact) {}

  BiCGSTAB(ParameterList &parlist, bool useInexact = false)
    : Krylov<Real>(parlist), isInitialized_(false), useInexact_(useInexact) {}

  Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M,
            int &iter, int &flag ) {
    if ( !isInitialized_ ) {
      r_ = b.clone(); r1_ = b.clone(); p_ = b.clone();
      v_ = b.clone(); s_  = b.clone(); t_ = b.clone();
      h_ = x.clone(); y_  = x.clone(); z_ = x.clone();
      isInitialized_ = true;
    }

    Real rho(1), rho1(1), alpha(1), beta(0), omega(1);
    Real rnorm = b.norm();
    Real itol = std::sqrt(ROL_EPSILON<Real>());
    const Real rtol = std::min(Krylov<Real>::getAbsoluteTolerance(),Krylov<Real>::getRelativeTolerance()*rnorm);
    if (rnorm <= rtol) return rnorm;

    x.zero();
    v_->zero();
    p_->zero();
    r_->set(b);
    r1_->set(*r_);

    iter = 0;
    flag = 0;

    for (iter = 0; iter < (int)Krylov<Real>::getMaximumIteration(); iter++) {
      rho1 = r_->dot(*r1_);
      beta = (rho1/rho)*(alpha/omega);
      p_->axpy(-omega,*v_);
      p_->scale(beta);
      p_->plus(*r_);

      if ( useInexact_ )
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm);
      M.applyInverse(*y_, *p_, itol);
      A.apply(*v_, *y_, itol);
      
      alpha = rho1 / v_->dot(*r1_);
      h_->set(x);
      h_->axpy(alpha,*y_);
      s_->set(*r_);
      s_->axpy(-alpha,*v_);

      rnorm = s_->norm();
      if (rnorm <= rtol) {
        x.set(*h_);
        break;
      }

      if ( useInexact_ )
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm);
      M.applyInverse(*z_, *s_, itol);
      A.apply(*t_, *z_, itol);

      omega = t_->dot(*s_) / t_->dot(*t_);
      x.set(*h_);
      x.axpy(omega,*z_);
      r_->set(*s_);
      r_->axpy(-omega,*t_);

      rnorm = r_->norm();
      if (rnorm <= rtol) break;

      rho = rho1;
    }
    if (iter == (int)Krylov<Real>::getMaximumIteration()) {
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
