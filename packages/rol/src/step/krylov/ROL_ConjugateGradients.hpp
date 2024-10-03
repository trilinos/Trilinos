// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONJUGATEGRADIENTS_H
#define ROL_CONJUGATEGRADIENTS_H

/** \class ROL::ConjugateGradients
    \brief Provides definitions of the Conjugate Gradient solver.
*/

#include "ROL_Krylov.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class ConjugateGradients : public Krylov<Real> {

  bool isInitialized_;
  bool useInexact_;
  ROL::Ptr<Vector<Real> > r_;
  ROL::Ptr<Vector<Real> > v_;
  ROL::Ptr<Vector<Real> > p_;
  ROL::Ptr<Vector<Real> > Ap_;

public:
  ConjugateGradients(Real absTol = 1.e-4, Real relTol = 1.e-2, unsigned maxit = 100, bool useInexact = false)
    : Krylov<Real>(absTol,relTol,maxit), isInitialized_(false), useInexact_(useInexact) {}

  Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M,
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
    Real itol = std::sqrt(ROL_EPSILON<Real>());

    x.zero();
    r_->set(b);

    M.applyInverse(*v_, *r_, itol);
    p_->set(*v_);

    iter = 0;
    flag = 0;

    Real kappa(0), beta(0), alpha(0), tmp(0), zero(0);
    //Real gv    = v_->dot(r_->dual());
    Real gv    = v_->apply(*r_);

    for (iter = 0; iter < (int)Krylov<Real>::getMaximumIteration(); iter++) {
      if ( useInexact_ ) {
        itol = rtol/((Real)Krylov<Real>::getMaximumIteration() * rnorm);
      }
      A.apply(*Ap_, *p_, itol);

      //kappa = p_->dot(Ap_->dual());
      kappa = p_->apply(*Ap_);
      if ( kappa <= zero ) {
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

      itol = std::sqrt(ROL_EPSILON<Real>());
      M.applyInverse(*v_, *r_, itol);
      tmp  = gv;
      //gv   = v_->dot(r_->dual());
      gv   = v_->apply(*r_);
      beta = gv/tmp;
 
      p_->scale(beta);
      p_->plus(*v_);
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
