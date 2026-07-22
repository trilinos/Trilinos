// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LDFP_H
#define ROL_LDFP_H

/** \class ROL::lDFP
    \brief Provides definitions for limited-memory DFP operators.
*/

#include "ROL_Secant.hpp"

namespace ROL {

template<class Real>
class lDFP : public Secant<Real> {
private:
  using Secant<Real>::state_;
  using Secant<Real>::useDefaultScaling_;
  using Secant<Real>::Bscaling_;

public:
  lDFP(int M, bool useDefaultScaling = true, Real Bscaling = Real(1))
    : Secant<Real>(M,useDefaultScaling,Bscaling) {}

  // Apply lBFGS Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const {
    const Real one(1);

    // Apply initial Hessian approximation to v
    applyH0(Hv,v);

    std::vector<Ptr<Vector<Real>>> a(state_->current+1);
    std::vector<Ptr<Vector<Real>>> b(state_->current+1);
    Real bv(0), av(0), bs(0), as(0);
    for (int i = 0; i <= state_->current; i++) {
      b[i] = Hv.clone();
      b[i]->set(*(state_->iterDiff[i]));
      b[i]->scale(1.0/sqrt(state_->product[i]));
      //bv = b[i]->dot(v.dual());
      bv = b[i]->apply(v);
      Hv.axpy(bv,*b[i]);

      a[i] = Hv.clone();
      applyH0(*a[i],*(state_->gradDiff[i]));

      for (int j = 0; j < i; j++) {
        //bs = b[j]->dot((state_->gradDiff[i])->dual());
        bs = b[j]->apply(*(state_->gradDiff[i]));
        a[i]->axpy(bs,*b[j]);
        //as = a[j]->dot((state_->gradDiff[i])->dual());
        as = a[j]->apply(*(state_->gradDiff[i]));
        a[i]->axpy(-as,*a[j]);
      }
      //as = a[i]->dot((state_->gradDiff[i])->dual());
      as = a[i]->apply(*(state_->gradDiff[i]));
      a[i]->scale(one/sqrt(as));
      //av = a[i]->dot(v.dual());
      av = a[i]->apply(v);
      Hv.axpy(-av,*a[i]);
    }
  }

  // Apply Initial Secant Approximate Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v ) const {
    Hv.set(v.dual());
    if (useDefaultScaling_) {
      if (state_->iter != 0 && state_->current != -1) {
        Real ss = state_->iterDiff[state_->current]->dot(*(state_->iterDiff[state_->current]));
        Hv.scale(state_->product[state_->current]/ss);
      }
    }
    else {
      Hv.scale(static_cast<Real>(1)/Bscaling_);
    }
  }

  // Apply lBFGS Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const {
    const Real zero(0);

    Bv.set(v.dual());
    std::vector<Real> alpha(state_->current+1,zero);
    for (int i = state_->current; i>=0; i--) {
      alpha[i]  = state_->gradDiff[i]->dot(Bv);
      alpha[i] /= state_->product[i];
      Bv.axpy(-alpha[i],(state_->iterDiff[i])->dual());
    }

    // Apply initial inverse Hessian approximation to v
    Ptr<Vector<Real>> tmp = Bv.clone();
    applyB0(*tmp,Bv.dual());
    Bv.set(*tmp);

    Real beta(0);
    for (int i = 0; i <= state_->current; i++) {
      //beta  = state_->iterDiff[i]->dot(Bv.dual());
      beta  = state_->iterDiff[i]->apply(Bv);
      beta /= state_->product[i];
      Bv.axpy((alpha[i]-beta),*(state_->gradDiff[i]));
    }
  }

  // Apply Initial Secant Approximate Hessian
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v ) const {
    Bv.set(v.dual());
    if (useDefaultScaling_) {
      if (state_->iter != 0 && state_->current != -1) {
        Real ss = state_->iterDiff[state_->current]->dot(*(state_->iterDiff[state_->current]));
        Bv.scale(ss/state_->product[state_->current]);
      }
    }
    else {
      Bv.scale(Bscaling_);
    }
  }
};

}

#endif
