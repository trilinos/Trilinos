// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LSR1_H
#define ROL_LSR1_H

/** \class ROL::lSR1
    \brief Provides definitions for limited-memory SR1 operators.
*/

#include "ROL_Secant.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class lSR1 : public Secant<Real> {
private:

  //mutable bool updateIterate_;
  bool isInitialized_;
  mutable bool H0called_, B0called_;
  Ptr<Vector<Real>> Bs_, Hy_, prim_, dual_;

  using Secant<Real>::state_;
  using Secant<Real>::y_;
  using Secant<Real>::useDefaultScaling_;
  using Secant<Real>::Bscaling_;

public:
  lSR1(int M, bool useDefaultScaling = true, Real Bscaling = Real(1), ESecantMode mode = SECANTMODE_BOTH)
    : Secant<Real>(M,useDefaultScaling,Bscaling,mode), isInitialized_(false),
      H0called_(false), B0called_(false) {
    if (useDefaultScaling_) Bscaling_ = static_cast<Real>(1);
    //updateIterate_ = true;
  }

  // Update Secant Approximation
  void updateStorage( const Vector<Real> &x,  const Vector<Real> &grad,
                      const Vector<Real> &gp, const Vector<Real> &s,
                      const Real snorm,       const int iter ) {
    const Real one(1), tol(std::sqrt(ROL_EPSILON<Real>()));
    if ( !isInitialized_ ) {
      state_->iterate = x.clone();
      y_              = grad.clone();
      if (state_->mode == SECANTMODE_FORWARD) {
        Bs_ = grad.clone(); dual_ = grad.clone();
      }
      else if (state_->mode == SECANTMODE_INVERSE) {
        Hy_ = x.clone();    prim_ = x.clone();
      }
      else {
        Bs_ = grad.clone(); dual_ = grad.clone();
        Hy_ = x.clone();    prim_ = x.clone();
      }
      isInitialized_ = true;
    }

    // Update iterate
    state_->iter = iter;
    state_->iterate->set(x);

    // Compute gradient difference
    y_->set(grad);
    y_->axpy(-one,gp);

    Real dotF(ROL_INF<Real>()), tolF(0), dotI(ROL_INF<Real>()), tolI(0);
    if (state_->mode == SECANTMODE_FORWARD || state_->mode == SECANTMODE_BOTH) {
      // Compute y - Bs and <s, y - Bs>
      applyB(*Bs_,s);
      Bs_->scale(-one);
      Bs_->plus(*y_);
      //dotF = s.dot(Bs_->dual());
      dotF = s.apply(*Bs_);
      tolF = tol*snorm*Bs_->norm();
    }
    if (state_->mode == SECANTMODE_INVERSE || state_->mode == SECANTMODE_BOTH) {
      // Compute s - Hy and <y, s - Hy>
      applyH(*Hy_,*y_);
      Hy_->scale(-one);
      Hy_->plus(s);
      //dotI = y_->dot(Hy_->dual());
      dotI = y_->apply(*Hy_);
      tolI = tol*y_->norm()*Hy_->norm();
    }
    if (std::abs(dotF) > tolF && std::abs(dotI) > tolI) {
      if (state_->current < state_->storage-1) {
        state_->current++;
        if (state_->mode == SECANTMODE_INVERSE || state_->mode == SECANTMODE_BOTH) {
          state_->iterDiff.push_back(x.clone());            // Create new memory
        }
        if (state_->mode == SECANTMODE_FORWARD || state_->mode == SECANTMODE_BOTH) {
          state_->gradDiff.push_back(grad.clone());         // Create new memory
        }
      }
      else {
        if (state_->mode == SECANTMODE_INVERSE || state_->mode == SECANTMODE_BOTH) {
          state_->iterDiff.push_back(state_->iterDiff[0]);  // Move first element to the last
          state_->iterDiff.erase(state_->iterDiff.begin()); // Remove first element of s list 
          state_->product2.erase(state_->product2.begin()); // Remove first element of rho list
        }
        if (state_->mode == SECANTMODE_FORWARD || state_->mode == SECANTMODE_BOTH) {
          state_->gradDiff.push_back(state_->gradDiff[0]);  // Move first element to the last
          state_->gradDiff.erase(state_->gradDiff.begin()); // Remove first element of y list
          state_->product.erase(state_->product.begin());   // Remove first element of rho list
        }
      }
      if (state_->mode == SECANTMODE_INVERSE || state_->mode == SECANTMODE_BOTH) {
        state_->iterDiff[state_->current]->set(*Hy_);       // s_k - H_k y_k
        state_->product2.push_back(dotI);                   // (s_k - H_k y_k)' y_k
      }
      if (state_->mode == SECANTMODE_FORWARD || state_->mode == SECANTMODE_BOTH) {
        state_->gradDiff[state_->current]->set(*Bs_);       // y_k - B_k s_k
        state_->product.push_back(dotF);                    // (y_k - B_k s_k)' s_k  
      }
      //if (useDefaultScaling_) Bscaling_ = s.dot(y_->dual())/(snorm*snorm);
      if (useDefaultScaling_) Bscaling_ = s.apply(*y_)/(snorm*snorm);
    }
    /*
    const Real one(1);
    if ( !isInitialized_ ) {
      state_->iterate = x.clone();
      y_             = grad.clone();
      isInitialized_ = true;
    }

    state_->iterate->set(x);
    state_->iter = iter;
    y_->set(grad);
    y_->axpy(-one,gp);

    if (updateIterate_ || state_->current == -1) {
      Real sy = s.dot(y_->dual());
      if (state_->current < state_->storage-1) {
        state_->current++;                                 // Increment Storage
        state_->iterDiff.push_back(s.clone());            // Create new memory
        state_->gradDiff.push_back(grad.clone());         // Create new memory
      }
      else {
        state_->iterDiff.push_back(state_->iterDiff[0]);  // Move first element to the last
        state_->gradDiff.push_back(state_->gradDiff[0]);  // Move first element to the last
        state_->iterDiff.erase(state_->iterDiff.begin());   // Remove first element of s list 
        state_->gradDiff.erase(state_->gradDiff.begin());   // Remove first element of y list
        state_->product.erase(state_->product.begin());     // Remove first element of rho list
      }
      state_->iterDiff[state_->current]->set(s);            // s=x_{k+1}-x_k
      state_->gradDiff[state_->current]->set(*y_);          // y=g_{k+1}-g_k
      state_->product.push_back(sy);                       // ys=1/rho  
    }
    updateIterate_ = true;
    */
  }

  // Apply Initial Secant Approximate Inverse Hessian
  virtual void applyH0( Vector<Real> &Hv, const Vector<Real> &v ) const {
    if (state_->current > -1) {
      prim_->set(v.dual());
      Hv.set(*prim_);
      H0called_ = true;
    }
    else {
      Hv.set(v.dual());
    }
    Hv.scale(static_cast<Real>(1)/Bscaling_);
  }

  // Apply lSR1 Approximate Inverse Hessian
  void applyH( Vector<Real> &Hv, const Vector<Real> &v ) const {
    if (state_->mode == SECANTMODE_INVERSE || state_->mode == SECANTMODE_BOTH) {
      // Apply initial Hessian approximation to v
      H0called_ = false;
      applyH0(Hv,v);
      // Apply rank one updates
      if (state_->current > -1) {
        Real prod(0);
        if (!H0called_) prim_->set(v.dual());
        for (int i = 0; i <= state_->current; ++i) {
          prod = state_->iterDiff[i]->dot(*prim_);
          Hv.axpy(prod/state_->product2[i],*state_->iterDiff[i]);
        }
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::lSR1::applyH : Not supported in forward mode!");
    }
    /*
    std::vector<Ptr<Vector<Real>>> a(state_->current+1);
    std::vector<Ptr<Vector<Real>>> b(state_->current+1);
    Real byi(0), byj(0), bv(0), normbi(0), normyi(0), one(1);
    for (int i = 0; i <= state_->current; i++) {
      // Compute Hy
      a[i] = Hv.clone();
      applyH0(*(a[i]),*(state_->gradDiff[i]));
      for (int j = 0; j < i; j++) {
        byj = b[j]->dot((state_->gradDiff[j])->dual());
        byi = b[j]->dot((state_->gradDiff[i])->dual());
        a[i]->axpy(byi/byj,*(b[j]));
      }
      // Compute s - Hy
      b[i] = Hv.clone();
      b[i]->set(*(state_->iterDiff[i]));
      b[i]->axpy(-one,*(a[i]));

      // Compute Hv
      byi    = b[i]->dot((state_->gradDiff[i])->dual());
      normbi = b[i]->norm();
      normyi = (state_->gradDiff[i])->norm();
      if ( i == state_->current && std::abs(byi) < sqrt(ROL_EPSILON<Real>())*normbi*normyi ) {
        updateIterate_ = false;
      }
      else {
        updateIterate_ = true;
        bv  = b[i]->dot(v.dual());
        Hv.axpy(bv/byi,*(b[i]));
      }
    }
    */
  }

  // Apply Initial Secant Approximate Hessian
  virtual void applyB0( Vector<Real> &Bv, const Vector<Real> &v ) const {
    if (state_->current > -1) {
      dual_->set(v.dual());
      Bv.set(*dual_);
      B0called_ = true;
    }
    else {
      Bv.set(v.dual());
    }
    Bv.scale(Bscaling_);
  }

  // Apply lSR1 Approximate Hessian
  void applyB( Vector<Real> &Bv, const Vector<Real> &v ) const {
    if (state_->mode == SECANTMODE_FORWARD || state_->mode == SECANTMODE_BOTH) {
      // Apply initial Hessian approximation to v
      B0called_ = false;
      applyB0(Bv,v);
      // Apply rank one updates
      if (state_->current > -1) {
        Real prod(0);
        if (!B0called_) dual_->set(v.dual());
        for (int i = 0; i <= state_->current; ++i) {
          prod = state_->gradDiff[i]->dot(*dual_);
          Bv.axpy(prod/state_->product[i],*state_->gradDiff[i]);
        }
      }
    }
    else {
      throw Exception::NotImplemented(">>> ROL::lSR1::applyB : Not supported in inverse mode!");
    }
    /*
    std::vector<Ptr<Vector<Real>>> a(state_->current+1);
    std::vector<Ptr<Vector<Real>>> b(state_->current+1);
    Real bsi(0), bsj(0), bv(0), normbi(0), normsi(0), one(1);
    for (int i = 0; i <= state_->current; i++) {
      // Compute Hy
      a[i] = Bv.clone();
      applyB0(*(a[i]),*(state_->iterDiff[i]));
      for (int j = 0; j < i; j++) {
        bsj = (state_->iterDiff[j])->dot(b[j]->dual());
        bsi = (state_->iterDiff[i])->dot(b[j]->dual());
        a[i]->axpy(bsi/bsj,*(b[j]));
      }
      // Compute s - Hy
      b[i] = Bv.clone();
      b[i]->set(*(state_->gradDiff[i]));
      b[i]->axpy(-one,*(a[i]));

      // Compute Hv
      bsi    = (state_->iterDiff[i])->dot(b[i]->dual());
      normbi = b[i]->norm();
      normsi = (state_->iterDiff[i])->norm();
      if ( i == state_->current && std::abs(bsi) < sqrt(ROL_EPSILON<Real>())*normbi*normsi ) {
        updateIterate_ = false;
      }
      else {
        updateIterate_ = true;
        bv  = b[i]->dot(v.dual());
        Bv.axpy(bv/bsi,*(b[i]));
      }
    }
    */
  }
};

}

#endif
