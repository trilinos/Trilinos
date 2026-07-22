// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PD_MEANSEMIDEVIATION_HPP
#define ROL_PD_MEANSEMIDEVIATION_HPP

#include "ROL_PD_RandVarFunctional.hpp"

namespace ROL {

template<class Real>
class PD_MeanSemiDeviation : public PD_RandVarFunctional<Real> {
private:
  Real coeff_;

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
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ < zero), std::invalid_argument,
      ">>> ERROR (ROL::PD_MeanSemiDeviation): Element of coefficient array out of range!");
    initializeStorage();
  }

public:
  PD_MeanSemiDeviation(const Real coeff)
    : PD_RandVarFunctional<Real>(), coeff_(coeff) {
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
    Real val = computeValue(obj,x,tol);
    val_ += weight_ * val;
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    // Compute expected value
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    // Compute deviation
    Real diff(0), pf0(0), dev(0), weight(0), lam(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      diff  -= ev;
      setValue(diff, sampler.getMyPoint(i));
      getMultiplier(lam,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      pf0   += weight * ppf(diff, lam, getPenaltyParameter(), 0);
    }
    sampler.sumAll(&pf0,&dev,1);
    dev *= coeff_;
    // Return mean plus deviation
    return ev + dev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_ += weight_ * val;
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    // Compute expected value
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    // Compute deviation
    hv_->zero(); dualVector_->zero();
    Real diff(0), pf(0), pf1(0), dev(0), one(1), weight(0), lam(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      diff  -= ev;
      getMultiplier(lam,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      pf1    = weight * ppf(diff, lam, getPenaltyParameter(), 1);
      pf    += pf1;
      gradients_->get(*hv_, sampler.getMyPoint(i));
      dualVector_->axpy(coeff_ * pf1, *hv_);
    }
    sampler.sumAll(&pf,&dev,1);
    g_->scale(one - coeff_ * dev);
    g_->plus(*dualVector_);
    sampler.sumAll(*g_,g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * val;
    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    gv_     += weight_ * gv;
    g_->axpy(weight_, *dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_, *dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real one(1);
    // Compute expected value
    std::vector<Real> mval(2), gval(2);
    mval[0] = val_; mval[1] = gv_;
    sampler.sumAll(&mval[0],&gval[0],2);
    Real ev = gval[0], egv = gval[1];
    // Compute deviation
    std::vector<Real> mvec(3), gvec(3);
    Real diff(0), gv(0), pf1(0), pf2(0), weight(0), lam(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      gradvecs_->get(gv,sampler.getMyPoint(i));
      getMultiplier(lam,sampler.getMyPoint(i));
      weight   = sampler.getMyWeight(i);
      diff    -= ev;
      pf1      = ppf(diff, lam, getPenaltyParameter(), 1);
      pf2      = ppf(diff, lam, getPenaltyParameter(), 2);
      mvec[0] += weight * pf1;
      mvec[1] += weight * pf2;
      mvec[2] += weight * pf2 * gv;
    }
    sampler.sumAll(&mvec[0],&gvec[0],3);
    Real c1 = one - coeff_ * gvec[0];
    Real c2 = coeff_ * (gvec[1]*egv - gvec[2]);
    hv_->scale(c1);
    hv_->axpy(c2, *g_); 
    sampler.sumAll(*hv_,hv);

    dualVector_->zero(); hv_->zero(); g_->zero();
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      gradients_->get(*g_,sampler.getMyPoint(i));
      gradvecs_->get(gv,sampler.getMyPoint(i));
      hessvecs_->get(*dualVector_,sampler.getMyPoint(i));
      getMultiplier(lam,sampler.getMyPoint(i));
      weight   = sampler.getMyWeight(i);
      diff    -= ev;
      pf1      = ppf(diff, lam, getPenaltyParameter(), 1);
      pf2      = ppf(diff, lam, getPenaltyParameter(), 2);
      hv_->axpy(weight * coeff_ * pf2 * (gv-egv), *g_);
      hv_->axpy(weight * coeff_ *pf1, *dualVector_);
    }
    sampler.sumAll(*hv_, *dualVector_);
    hv.plus(*dualVector_);
  }
};

}

#endif
