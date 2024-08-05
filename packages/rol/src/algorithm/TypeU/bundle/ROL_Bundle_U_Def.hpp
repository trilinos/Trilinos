// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_U_DEF_H
#define ROL_BUNDLE_U_DEF_H

#include "ROL_Types.hpp"

namespace ROL {

template<typename Real>
void Bundle_U<Real>::remove(const std::vector<unsigned> &ind) {
  Real zero(0);
  for (unsigned j = ind.back()+1; j < size_; ++j) {
    (subgradients_[j-1])->set(*(subgradients_[j]));
    linearizationErrors_[j-1] = linearizationErrors_[j];
    distanceMeasures_[j-1]    = distanceMeasures_[j];
    dualVariables_[j-1]       = dualVariables_[j];
  }
  (subgradients_[size_-1])->zero();
  linearizationErrors_[size_-1] = ROL_OVERFLOW<Real>();
  distanceMeasures_[size_-1]    = ROL_OVERFLOW<Real>();
  dualVariables_[size_-1]       = zero;
  for (unsigned i = ind.size()-1; i > 0; --i) {
    for (unsigned j = ind[i-1]+1; j < size_; ++j) {
      (subgradients_[j-1])->set(*(subgradients_[j]));
      linearizationErrors_[j-1] = linearizationErrors_[j];
      distanceMeasures_[j-1]    = distanceMeasures_[j];
      dualVariables_[j-1]       = dualVariables_[j];
    }
  }
  size_ -= ind.size();
}

template<typename Real>
void Bundle_U<Real>::add(const Vector<Real> &g, const Real le, const Real dm) {
  Real zero(0);
  (subgradients_[size_])->set(g);
  linearizationErrors_[size_] = le;
  distanceMeasures_[size_]    = dm;
  dualVariables_[size_]       = zero;
  size_++;
}
 
template<typename Real>
Bundle_U<Real>::Bundle_U(const unsigned maxSize,
                         const Real coeff,
                         const Real omega,
                         const unsigned remSize) 
  : size_(0), maxSize_(maxSize), remSize_(remSize),
    coeff_(coeff), omega_(omega), isInitialized_(false) {
  Real zero(0);
  remSize_ = ((remSize_ < 2) ? 2 : ((remSize_ > maxSize_-1) ? maxSize_-1 : remSize_));
  coeff_ = std::max(static_cast<Real>(0),coeff_);
  omega_ = std::max(static_cast<Real>(1),omega_);
  subgradients_.clear();
  subgradients_.assign(maxSize,nullPtr);
  linearizationErrors_.clear();
  linearizationErrors_.assign(maxSize_,ROL_OVERFLOW<Real>());
  distanceMeasures_.clear();
  distanceMeasures_.assign(maxSize_,ROL_OVERFLOW<Real>());
  dualVariables_.clear();
  dualVariables_.assign(maxSize_,zero);
}

template<typename Real>
void Bundle_U<Real>::initialize(const Vector<Real> &g) {
  if ( !isInitialized_ ) {
    Real zero(0), one(1);
    for (unsigned i = 0; i < maxSize_; ++i) {
      subgradients_[i] = g.clone();
    }
    (subgradients_[0])->set(g);
    linearizationErrors_[0] = zero;
    distanceMeasures_[0]    = zero;
    dualVariables_[0]       = one;
    size_++;
    isInitialized_ = true;
    tG_ = g.clone();
    yG_ = g.clone();
    eG_ = g.clone();
    gx_ = g.clone();
    ge_ = g.clone();
  }
}

template<typename Real>
const Real Bundle_U<Real>::linearizationError(const unsigned i) const {
  return linearizationErrors_[i];
}

template<typename Real>
const Real Bundle_U<Real>::distanceMeasure(const unsigned i) const {
  return distanceMeasures_[i];
}

template<typename Real>
const Vector<Real> & Bundle_U<Real>::subgradient(const unsigned i) const {
  return *(subgradients_[i]);
}
  
template<typename Real>
const Real Bundle_U<Real>::getDualVariable(const unsigned i) const {
  return dualVariables_[i];
}
  
template<typename Real>
void Bundle_U<Real>::setDualVariable(const unsigned i, const Real val) {
  dualVariables_[i] = val;
}

template<typename Real>
void Bundle_U<Real>::resetDualVariables(void) {
  const Real zero(0);
  dualVariables_.assign(size_,zero);
}

template<typename Real>
const Real Bundle_U<Real>::computeAlpha(const Real dm, const Real le) const {
  Real alpha = le;
  if ( coeff_ > ROL_EPSILON<Real>() ) {
    alpha = std::max(coeff_*std::pow(dm,omega_),le);
  }
  return alpha;
}

template<typename Real>
const Real Bundle_U<Real>::alpha(const unsigned i) const {
  return computeAlpha(distanceMeasures_[i],linearizationErrors_[i]);
}

template<typename Real>
unsigned Bundle_U<Real>::size(void) const {
  return size_;
}

template<typename Real>
void Bundle_U<Real>::aggregate(Vector<Real> &aggSubGrad, Real &aggLinErr, Real &aggDistMeas) const {
  Real zero(0), one(1);
  aggSubGrad.zero(); aggLinErr = zero; aggDistMeas = zero; eG_->zero();
  Real eLE(0), eDM(0), yLE(0), yDM(0), tLE(0), tDM(0);
  for (unsigned i = 0; i < size_; ++i) {
    // Compute aggregate subgradient using Kahan's compensated sum
    //aggSubGrad.axpy(dualVariables_[i],*subgradients_[i]);
    yG_->set(*subgradients_[i]); yG_->scale(dualVariables_[i]); yG_->axpy(-one,*eG_);
    tG_->set(aggSubGrad); tG_->plus(*yG_);
    eG_->set(*tG_); eG_->axpy(-one,aggSubGrad); eG_->axpy(-one,*yG_);
    aggSubGrad.set(*tG_);
    // Compute aggregate linearization error using Kahan's compensated sum
    //aggLinErr += dualVariables_[i]*linearizationErrors_[i];
    yLE = dualVariables_[i]*linearizationErrors_[i] - eLE;
    tLE = aggLinErr + yLE;
    eLE = (tLE - aggLinErr) - yLE;
    aggLinErr = tLE;
    // Compute aggregate distance measure using Kahan's compensated sum
    //aggDistMeas += dualVariables_[i]*distanceMeasures_[i];
    yDM = dualVariables_[i]*distanceMeasures_[i] - eDM;
    tDM = aggDistMeas + yDM;
    eDM = (tDM - aggDistMeas) - yDM;
    aggDistMeas = tDM;
  }
}

template<typename Real>
void Bundle_U<Real>::reset(const Vector<Real> &g, const Real le, const Real dm) {
  if (size_ == maxSize_) {
    // Find indices to remove
    unsigned loc = size_, cnt = 0;
    std::vector<unsigned> ind(remSize_,0);
    for (unsigned i = size_; i > 0; --i) {
      if ( std::abs(linearizationErrors_[i-1]) < ROL_EPSILON<Real>() ) {
        loc = i-1;
        break;
      }
    }
    for (unsigned i = 0; i < size_; ++i) {
      if ( i != loc ) {
        ind[cnt] = i;
        cnt++;
      }
      if (cnt == remSize_) {
        break;
      }
    }
    // Remove indices
    remove(ind);
    // Add aggregate subgradient
    add(g,le,dm);
  }
}

template<typename Real>
void Bundle_U<Real>::update(const bool flag, const Real linErr, const Real distMeas,
                            const Vector<Real> &g, const Vector<Real> &s) {
  Real zero(0);
  if ( flag ) {
    // Serious step taken: Update linearlization errors and distance measures
    for (unsigned i = 0; i < size_; ++i) {
      //linearizationErrors_[i] += linErr - subgradients_[i]->dot(s.dual());
      linearizationErrors_[i] += linErr - subgradients_[i]->apply(s);
      distanceMeasures_[i]    += distMeas;
    }
    linearizationErrors_[size_] = zero;
    distanceMeasures_[size_]    = zero;
  }
  else {
    // Null step taken: 
    linearizationErrors_[size_] = linErr;
    distanceMeasures_[size_]    = distMeas;
  }
  // Update (sub)gradient bundle
  (subgradients_[size_])->set(g);
  // Update dual variables
  dualVariables_[size_] = zero;
  // Update bundle size
  size_++;
}

template<typename Real>
const Real Bundle_U<Real>::GiGj(const unsigned i, const unsigned j) const {
  return subgradient(i).dot(subgradient(j));
}

template<typename Real>
const Real Bundle_U<Real>::dotGi(const unsigned i, const Vector<Real> &x) const {
  return x.dot(subgradient(i));
}

template<typename Real>
void Bundle_U<Real>::addGi(const unsigned i, const Real a, Vector<Real> &x) const {
  x.axpy(a,subgradient(i));
}

template<typename Real>
Real Bundle_U<Real>::evaluateObjective(std::vector<Real> &g, const std::vector<Real> &x, const Real t) const {
  Real one(1), half(0.5);
  gx_->zero(); eG_->zero();
  for (unsigned i = 0; i < size(); ++i) {
    // Compute Gx using Kahan's compensated sum
    //gx_->axpy(x[i],*subgradients_[i]);
    yG_->set(subgradient(i)); yG_->scale(x[i]); yG_->axpy(-one,*eG_);
    tG_->set(*gx_); tG_->plus(*yG_);
    eG_->set(*tG_); eG_->axpy(-one,*gx_); eG_->axpy(-one,*yG_);
    gx_->set(*tG_);
  }
  Real Hx(0), val(0), err(0), tmp(0), y(0);
  for (unsigned i = 0; i < size(); ++i) {
    // Compute < g_i, Gx >
    Hx   = dotGi(i,*gx_);
    // Add to the objective function value using Kahan's compensated sum
    //val += x[i]*(half*Hx + alpha(i)/t);
    y    = x[i]*(half*Hx + alpha(i)/t) - err;
    tmp  = val + y;
    err  = (tmp - val) - y;
    val  = tmp;
    // Add gradient component
    g[i] = Hx + alpha(i)/t;
  }
  return val;
}

template<typename Real>
unsigned Bundle_U<Real>::solveDual_dim1(const Real t, const unsigned maxit, const Real tol) {
  setDualVariable(0,static_cast<Real>(1));
  //std::cout << "dim = " << Bundle<Real>::size() << "  iter = " << 0 << "  CONVERGED!\n";
  return 0;
}

template<typename Real>
unsigned Bundle_U<Real>::solveDual_dim2(const Real t, const unsigned maxit, const Real tol) {
  Real diffg  = gx_->dot(*gx_), zero(0), one(1), half(0.5);
  gx_->set(subgradient(0)); addGi(1,-one,*gx_);
  if ( std::abs(diffg) > ROL_EPSILON<Real>() ) {
    Real diffa  = (alpha(0)-alpha(1))/t;
    Real gdiffg = dotGi(1,*gx_);
    setDualVariable(0,std::min(one,std::max(zero,-(gdiffg+diffa)/diffg)));
    setDualVariable(1,one-getDualVariable(0));
  }
  else {
    if ( std::abs(alpha(0)-alpha(1)) > ROL_EPSILON<Real>() ) {
      if ( alpha(0) < alpha(1) ) {
        setDualVariable(0,one); setDualVariable(1,zero);
      }
      else if ( alpha(0) > alpha(1) ) {
        setDualVariable(0,zero); setDualVariable(1,one);
      }
    }
    else {
      setDualVariable(0,half); setDualVariable(1,half);
    }
  }
  //std::cout << "dim = " << Bundle<Real>::size() << "  iter = " << 0 << "  CONVERGED!\n";
  return 0;
}

} // namespace ROL

#endif
