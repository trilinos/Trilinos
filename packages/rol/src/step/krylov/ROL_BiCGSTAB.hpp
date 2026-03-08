// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
