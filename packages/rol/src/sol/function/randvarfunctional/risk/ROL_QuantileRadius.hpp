// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUANTILERADIUSQUADRANGLE_HPP
#define ROL_QUANTILERADIUSQUADRANGLE_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PlusFunction.hpp"

#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class QuantileRadius : public RandVarFunctional<Real> {
private:
  Ptr<PlusFunction<Real> > plusFunction_;
  Real prob_;
  Real coeff_;
  std::vector<Real> vec_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  void initializeQR(void) {
    Real zero(0);
    // Initialize temporary storage
    vec_.clear();  vec_.resize(2,zero);
  }

  void checkInputs(void) {
    Real zero(0), one(1);
    // Check inputs
    ROL_TEST_FOR_EXCEPTION((prob_>one || prob_<zero), std::invalid_argument,
      ">>> ERROR (ROL::QuantileRadius): Confidence level out of range!");
    ROL_TEST_FOR_EXCEPTION((coeff_<zero), std::invalid_argument,
      ">>> ERROR (ROL::QuantileRadius): Coefficient is negative!");
     initializeQR();
  }

public:

  QuantileRadius( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Quantile Radius");
    // Grab probability and coefficient arrays
    prob_  = list.get<Real>("Confidence Level");
    coeff_ = list.get<Real>("Coefficient");
    // Build (approximate) plus function
    plusFunction_ = makePtr<PlusFunction<Real>>(list);
    checkInputs();
  }

  QuantileRadius(const Real prob, const Real coeff,
                 const Ptr<PlusFunction<Real> > &pf)
    : RandVarFunctional<Real>(), plusFunction_(pf), prob_(prob), coeff_(coeff) {
    checkInputs();
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    vec_.assign(2,static_cast<Real>(0));
  }

  Real computeStatistic(const Ptr<const std::vector<Real>> &xstat) const override {
    Real stat(0), half(0.5);
    if (xstat != nullPtr) {
      stat = half*((*xstat)[0] + (*xstat)[1]);
    }
    return stat;
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    const Real half(0.5), one(1);
    Real val = computeValue(obj,x,tol);
    Real pf1 = plusFunction_->evaluate(val-xstat[0],0);
    Real pf2 = plusFunction_->evaluate(-val-xstat[1],0);
    RandVarFunctional<Real>::val_ += weight_*(val + half*coeff_/(one-prob_)*(pf1 + pf2));
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    const Real half(0.5);
    Real cvar(0);
    sampler.sumAll(&val_,&cvar,1);
    cvar += half*coeff_*(xstat[0] + xstat[1]);
    return cvar;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real half(0.5), one(1);
    Real val = computeValue(obj,x,tol);
    Real pf1 = plusFunction_->evaluate(val-xstat[0],1);
    Real pf2 = plusFunction_->evaluate(-val-xstat[1],1);
    Real c   = half*weight_*coeff_/(one-prob_);
    vec_[0] -= c*pf1;
    vec_[1] -= c*pf2;
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_ + c * (pf1 - pf2),*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real half(0.5);
    sampler.sumAll(&vec_[0],&gstat[0],2);
    sampler.sumAll(*g_,g);
    gstat[0] += half*coeff_;
    gstat[1] += half*coeff_;
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real half(0.5), one(1);
    Real val = computeValue(obj,x,tol);
    Real pf11 = plusFunction_->evaluate(val-xstat[0],1);
    Real pf12 = plusFunction_->evaluate(val-xstat[0],2);
    Real pf21 = plusFunction_->evaluate(-val-xstat[1],1);
    Real pf22 = plusFunction_->evaluate(-val-xstat[1],2);
    Real c    = half*weight_*coeff_/(one-prob_);
    Real gv   = computeGradVec(*dualVector_,obj,v,x,tol);
    vec_[0]  -= c*pf12*(gv-vstat[0]);
    vec_[1]  += c*pf22*(gv+vstat[1]);
    hv_->axpy(c*(pf12*(gv-vstat[0]) + pf22*(gv+vstat[1])),*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_ + c * (pf11 - pf21),*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(&vec_[0],&hvstat[0],2);
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
