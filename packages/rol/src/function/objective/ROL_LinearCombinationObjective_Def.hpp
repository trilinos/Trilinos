// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEARCOMBINATIONOBJECTIVE_DEF_H
#define ROL_LINEARCOMBINATIONOBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
LinearCombinationObjective<Real>::LinearCombinationObjective(const std::vector<Ptr<Objective<Real>>> &obj)
  : Objective<Real>(), obj_(obj),
    xdual_(nullPtr), initialized_(false) {
  size_ = obj_.size();
  weights_.clear(); weights_.assign(size_,static_cast<Real>(1));
}

template<typename Real>
LinearCombinationObjective<Real>::LinearCombinationObjective(const std::vector<Real> &weights,
                                                             const std::vector<Ptr<Objective<Real>>> &obj)
  : Objective<Real>(), obj_(obj), weights_(weights), size_(weights.size()),
    xdual_(nullPtr), initialized_(false) {}

template<typename Real>
void LinearCombinationObjective<Real>::update(const Vector<Real> &x, UpdateType type, int iter) {
  for (size_t i=0; i<size_; ++i) {
    obj_[i]->update(x,type,iter);
  }
}

template<typename Real>
void LinearCombinationObjective<Real>::update(const Vector<Real> &x, bool flag, int iter) {
  for (size_t i=0; i<size_; ++i) {
    obj_[i]->update(x,flag,iter);
  }
}

template<typename Real>
Real LinearCombinationObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  Real val(0);
  for (size_t i = 0; i < size_; i++) {
    val += weights_[i]*obj_[i]->value(x,tol);
  }
  return val;
}

template<typename Real>
void LinearCombinationObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  if (!initialized_) {
    xdual_ = g.clone();
    initialized_ = true;
  }
  g.zero();
  for (size_t i = 0; i < size_; i++) {
    obj_[i]->gradient(*xdual_,x,tol);
    g.axpy(weights_[i],*xdual_);
  }
}

template<typename Real>
void LinearCombinationObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  if (!initialized_) {
    xdual_ = hv.clone();
    initialized_ = true;
  }
  hv.zero();
  for (size_t i = 0; i < size_; i++) {
    obj_[i]->hessVec(*xdual_,v,x,tol);
    hv.axpy(weights_[i],*xdual_);
  }
}

template<typename Real>
void LinearCombinationObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  for (size_t i = 0; i < size_; ++i) {
    obj_[i]->setParameter(param);
  }
}

} // namespace ROL

#endif
