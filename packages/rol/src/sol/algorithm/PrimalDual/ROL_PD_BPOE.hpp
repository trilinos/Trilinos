// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PD_BPOE_HPP
#define ROL_PD_BPOE_HPP

#include "ROL_PD_RandVarFunctional.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class PD_BPOE : public PD_RandVarFunctional<Real> {
private:
  Real thresh_;

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
    initializeStorage();
  }

public:
  PD_BPOE(const Real thresh)
    : PD_RandVarFunctional<Real>(), thresh_(thresh) {
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
    Real arg = xstat[0] * (val - thresh_) + one;
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 0);
    val_  += weight_ * pf;
    setValue(arg, point_);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    return ev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real zero(0), one(1);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = xstat[0] * (val - thresh_) + one;
    Real  pf = ppf(arg, lam, getPenaltyParameter(), 1);
    if ( pf > zero ) {
      computeGradient(*dualVector_, obj, x, tol);
      val_ += weight_ * pf * (val - thresh_);
      g_->axpy(weight_ * pf * xstat[0], *dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    sampler.sumAll(*g_, g);
    gstat[0] = ev;
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
    Real arg = xstat[0] * (val - thresh_) + one;
    Real pf1 = ppf(arg, lam, getPenaltyParameter(), 1);
    Real pf2 = ppf(arg, lam, getPenaltyParameter(), 2);
    if ( pf1 > zero ) {
      Real gv = computeGradVec(*dualVector_, obj, v, x, tol);
      val_   += weight_ * pf1 * gv;
      hv_->axpy(weight_ * pf1 * vstat[0], *dualVector_);
      computeHessVec(*dualVector_, obj, v, x, tol);
      hv_->axpy(weight_ * pf1 * xstat[0], *dualVector_);
    }
    if ( pf2 > zero ) {
      Real gv = computeGradVec(*dualVector_, obj, v, x, tol);
      Real c1 = pf2 * (val - thresh_) * xstat[0];
      Real c2 = pf2 * (val - thresh_) * (val - thresh_);
      Real c3 = pf2 * xstat[0] * xstat[0];
      val_   += weight_ * (c1 * gv + c2 * vstat[0]);
      hv_->axpy(weight_ * (c3 * gv + c1 * vstat[0]), *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    sampler.sumAll(*hv_, hv);
    hvstat[0] = ev;
  }
};

}

#endif
