// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PH_PROBOBJECTIVE_H
#define PH_PROBOBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::PH_ProbObjective
    \brief Provides the interface for the progressive hedging probability objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_ProbObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Real threshold_;
  Real eps_;

  bool isValueComputed_;
  Real val_;

  bool isGradientInitialized_;
  bool isGradientComputed_;
  Ptr<Vector<Real>> g_;

  void getValue(const Vector<Real> &x, Real &tol) {
    if (!isValueComputed_) {
      val_ = obj_->value(x,tol);
      isValueComputed_ = true;
    }
  }

  void getGradient(const Vector<Real> &x, Real &tol) {
    if (!isGradientInitialized_) {
      g_ = x.dual().clone();
      isGradientInitialized_ = true;
    }
    if (!isGradientComputed_) {
      obj_->gradient(*g_,x,tol);
      isGradientComputed_ = true;
    }
  }

  Real smoothHeaviside(const Real x, const int deriv = 0) const {
    const Real one(1), two(2);
    Real val(0);
    if (deriv == 0) {
      Real ex = std::exp(-two*x/eps_);
      val = one/(one+ex);
    }
    else if (deriv == 1) {
      Real ex = std::exp(-two*x/eps_);
      val = (two/eps_)*ex/std::pow(one+ex,2);
    }
    else if (deriv == 2) {
      Real ex = std::exp(two*x/eps_);
      val = std::pow(two/eps_,2)*ex*(one-ex)/std::pow(one+ex,3);
    }
    return val;
  }

public:

  PH_ProbObjective(const Ptr<Objective<Real>> &obj,
                   ParameterList              &parlist)
    : obj_(obj),
      isValueComputed_(false),
      isGradientInitialized_(false),
      isGradientComputed_(false) {
    ParameterList &list = parlist.sublist("SOL").sublist("Probability").sublist("Smoothed POE");
    threshold_ = list.get<Real>("Threshold");
    eps_       = list.get<Real>("Smoothing Parameter");
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    isValueComputed_    = false;
    isGradientComputed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real prob = smoothHeaviside(val_-threshold_,0);
    if (std::abs(prob) > ROL_EPSILON<Real>()) {
      return prob;
    }
    return static_cast<Real>(0);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real prob = smoothHeaviside(val_-threshold_,1);
    if (std::abs(prob) > ROL_EPSILON<Real>()) {
      getGradient(x,tol);
      g.set(*g_); g.scale(prob);
    }
    else {
      g.zero();
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real prob1 = smoothHeaviside(val_-threshold_,1);
    Real prob2 = smoothHeaviside(val_-threshold_,2);
    if (std::abs(prob1) > ROL_EPSILON<Real>()) {
      obj_->hessVec(hv,v,x,tol);
      hv.scale(prob1);
    }
    else {
      hv.zero();
    }
    if (std::abs(prob2) > ROL_EPSILON<Real>()) {
      getGradient(x,tol);
      //Real gv    = v.dot(g_->dual());
      Real gv    = v.apply(*g_);
      hv.axpy(prob2*gv,*g_);
    }
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
