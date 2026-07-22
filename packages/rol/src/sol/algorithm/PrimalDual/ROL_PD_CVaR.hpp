// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PD_CVAR_HPP
#define ROL_PD_CVAR_HPP

#include "ROL_PD_RandVarFunctional.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class PD_CVaR : public PD_RandVarFunctional<Real> {
private:
  Real alpha_;
  Real beta_;

  Ptr<ScalarController<Real>> values_;
  Ptr<ScalarController<Real>> gradvecs_;
  Ptr<VectorController<Real>> gradients_;
  Ptr<VectorController<Real>> hessvecs_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  using PD_RandVarFunctional<Real>::setValue;
  using PD_RandVarFunctional<Real>::getMultiplier;
  using PD_RandVarFunctional<Real>::getPenaltyParameter;
  using PD_RandVarFunctional<Real>::ppf;

  void initializeStorage(void) {
    values_    = makePtr<ScalarController<Real>>();
    gradvecs_  = makePtr<ScalarController<Real>>();
    gradients_ = makePtr<VectorController<Real>>();
    hessvecs_  = makePtr<VectorController<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  void clear(void) {
    gradvecs_->reset();
    hessvecs_->reset();
  }

  void checkInputs(void) {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((alpha_ <= zero) || (alpha_ > one), std::invalid_argument,
      ">>> ERROR (ROL::PD_CVaR): Convex combination parameter alpha is out of range!");
    ROL_TEST_FOR_EXCEPTION((beta_ < zero) || (beta_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::PD_CVaR): Confidence parameter beta is out of range!");
    initializeStorage();
  }

public:
  PD_CVaR(const Real alpha, const Real beta)
    : PD_RandVarFunctional<Real>(), alpha_(alpha), beta_(beta) {
    checkInputs();
  }

  void setStorage(const Ptr<ScalarController<Real>> &value_storage,
                  const Ptr<VectorController<Real>> &gradient_storage) {
    values_    = value_storage;
    gradients_ = gradient_storage;
    PD_RandVarFunctional<Real>::setStorage(values_,gradients_);
  }

  void setHessVecStorage(const Ptr<ScalarController<Real>> &gradvec_storage,
                         const Ptr<VectorController<Real>> &hessvec_storage) {
    gradvecs_ = gradvec_storage;
    hessvecs_ = hessvec_storage;
    PD_RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }
 
  void initialize(const Vector<Real> &x) {
    PD_RandVarFunctional<Real>::initialize(x);
    clear();
  }
  
  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    const Real one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 0);
    val_    += weight_ * ((one-alpha_) * val + alpha_/(one-beta_) * pf);
    setValue(arg, point_);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    return ev + alpha_ * xstat[0];
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real zero(0), one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 1);
    val_    += weight_ * pf;
    Real c   = (one-alpha_) + alpha_/(one-beta_) * pf;
    if ( std::abs(c) > zero ) {
      computeGradient(*dualVector_, obj, x, tol);
      g_->axpy(weight_ * c, *dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real one(1);
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    ev *= -alpha_/(one-beta_);
    ev += alpha_;
    gstat[0] = ev;
    sampler.sumAll(*g_, g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real zero(0), one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real pf1 = ppf(arg, lam, getPenaltyParameter(), 1);
    Real pf2 = ppf(arg, lam, getPenaltyParameter(), 2);
    Real c(0);
    if ( std::abs(pf2) > zero ) {
      Real gv = computeGradVec(*dualVector_, obj, v, x, tol);
      val_   += weight_ * pf2 * (vstat[0] - gv);
      c       = pf2 * alpha_/(one-beta_) * (gv - vstat[0]);
      hv_->axpy(weight_ * c, *dualVector_);
    }
    c = (one-alpha_) + alpha_/(one-beta_) * pf1;
    if ( std::abs(c) > zero ) {
      computeHessVec(*dualVector_, obj, v, x, tol);
      hv_->axpy(weight_ * c, *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real one(1);
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    ev *= alpha_/(one-beta_);
    hvstat[0] = ev;
    sampler.sumAll(*hv_, hv);
  }
};

}

#endif
