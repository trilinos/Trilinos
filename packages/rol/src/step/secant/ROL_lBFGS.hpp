// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LBFGS_H
#define ROL_LBFGS_H

/** \class ROL::lBFGS
    \brief Provides definitions for limited-memory BFGS operators.
*/

#include "ROL_Secant.hpp"

namespace ROL {

template<class Real>
class lBFGS : public Secant<Real> {
private:
  using Secant<Real>::state_;

public:
  lBFGS(int M, bool useDefaultScaling = true, Real Bscaling = Real(1))
    : Secant<Real>(M,useDefaultScaling,Bscaling) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const {
    const Real zero(0);

    auto tmp = v.clone();
    tmp->set(v); 
    std::vector<Real> alpha(state_->current+1,zero);
    for (int i = state_->current; i>=0; i--) {
      alpha[i]  = state_->iterDiff[i]->apply(*tmp);
      alpha[i] /= state_->product[i];
      tmp->axpy(-alpha[i],*state_->gradDiff[i]);
    }

    // Apply initial inverse Hessian approximation to v
    Secant<Real>::applyH0(Hv,*tmp);

    Real beta(0);
    for (int i = 0; i <= state_->current; i++) {
      //beta  = Hv.dot((state_->gradDiff[i])->dual());
      beta  = Hv.apply(*state_->gradDiff[i]);
      beta /= state_->product[i];
      Hv.axpy((alpha[i]-beta),*(state_->iterDiff[i]));
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const {
    const Real one(1);

    // Apply initial Hessian approximation to v
    Secant<Real>::applyB0(Bv,v);

    std::vector<Ptr<Vector<Real>>> a(state_->current+1);
    std::vector<Ptr<Vector<Real>>> b(state_->current+1);
    Real bv(0), av(0), bs(0), as(0);
    for (int i = 0; i <= state_->current; i++) {
      b[i] = Bv.clone();
      b[i]->set(*(state_->gradDiff[i]));
      b[i]->scale(one/sqrt(state_->product[i]));
      //bv = v.dot(b[i]->dual());
      bv = v.apply(*b[i]);
      Bv.axpy(bv,*b[i]);

      a[i] = Bv.clone();
      Secant<Real>::applyB0(*a[i],*(state_->iterDiff[i]));

      for (int j = 0; j < i; j++) {
        //bs = (state_->iterDiff[i])->dot(b[j]->dual());
        bs = (state_->iterDiff[i])->apply(*b[j]);
        a[i]->axpy(bs,*b[j]);
        //as = (state_->iterDiff[i])->dot(a[j]->dual());
        as = (state_->iterDiff[i])->apply(*a[j]);
        a[i]->axpy(-as,*a[j]);
      }
      //as = (state_->iterDiff[i])->dot(a[i]->dual());
      as = (state_->iterDiff[i])->apply(*a[i]);
      a[i]->scale(one/sqrt(as));
      //av = v.dot(a[i]->dual());
      av = v.apply(*a[i]);
      Bv.axpy(-av,*a[i]);
    }
  }
};

}


#endif
