// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MIXEDQUANTILEQUADRANGLE_HPP
#define ROL_MIXEDQUANTILEQUADRANGLE_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PlusFunction.hpp"

#include "ROL_ParameterList.hpp"

/** @ingroup risk_group
    \class ROL::MixedCVaR
    \brief Provides an interface for a convex combination of
           conditional value-at-risks.

    The risk measure associated with the mixed-quantile quadrangle is defined
    as
    \f[
       \mathcal{R}(X) = \lambda_1 \mathrm{CVaR}_{\beta_1}(X)
         + \ldots + \lambda_n \mathrm{CVaR}_{\beta_n}(X)
    \f]
    where \f$0 \le \beta_1 \le \cdots \le \beta_n < 1\f$ and
    \f$0 \le \lambda_i\f$, \f$i=1,\ldots,n\f$, satisfies
    \f[
       \lambda_1 + \ldots + \lambda_n = 1.
    \f]
    Here, the conditional value-at-risk (CVaR) with confidence level
    \f$0\le \beta < 1\f$ is
    \f[
       \mathrm{CVaR}_\beta(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\beta} \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.  If the distribution of \f$X\f$ is
    continuous, then \f$\mathrm{CVaR}_{\beta}(X)\f$ is the conditional
    expectation of \f$X\f$ exceeding the \f$\beta\f$-quantile of \f$X\f$ and
    the optimal \f$t\f$ is the \f$\beta\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class MixedCVaR : public RandVarFunctional<Real> {
private:
  ROL::Ptr<PlusFunction<Real> > plusFunction_;
  std::vector<Real> prob_;
  std::vector<Real> coeff_;
  std::vector<Real> vec_;
  int size_;

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

  void initializeMCVAR(void) {
    size_ = prob_.size();
    vec_.clear();  vec_.resize(size_,static_cast<Real>(0));
  }

  void checkInputs(void) {
    int pSize = prob_.size(), cSize = coeff_.size();
    ROL_TEST_FOR_EXCEPTION((pSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MixedCVaR): Probability and coefficient arrays have different sizes!");
    Real sum(0), zero(0), one(1);
    for (int i = 0; i < pSize; i++) {
      ROL_TEST_FOR_EXCEPTION((prob_[i]>one || prob_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::MixedCVaR): Element of probability array out of range!");
      ROL_TEST_FOR_EXCEPTION((coeff_[i]>one || coeff_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::MixedCVaR): Element of coefficient array out of range!");
      sum += coeff_[i];
    }
    ROL_TEST_FOR_EXCEPTION((std::abs(sum-one) > std::sqrt(ROL_EPSILON<Real>())),std::invalid_argument,
      ">>> ERROR (ROL::MixedCVaR): Coefficients do not sum to one!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MixedCVaR): PlusFunction pointer is null!");
    initializeMCVAR();
  }

public:

  MixedCVaR( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mixed CVaR");
    // Grab probability and coefficient arrays
    prob_  = ROL::getArrayFromStringParameter<Real>(list,"Probability Array");
    coeff_ = ROL::getArrayFromStringParameter<Real>(list,"Coefficient Array");
    plusFunction_ = ROL::makePtr<PlusFunction<Real>>(list);
    // Check inputs
    checkInputs();
  }

  MixedCVaR(const std::vector<Real> &prob,
                          const std::vector<Real> &coeff,
                          const ROL::Ptr<PlusFunction<Real> > &pf )
    : RandVarFunctional<Real>(), plusFunction_(pf), prob_(prob), coeff_(coeff) {
    checkInputs();
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    vec_.assign(size_,static_cast<Real>(0));
  }

  Real computeStatistic(const Ptr<const std::vector<Real>> &xstat) const override {
    Real stat(0);
    if (xstat != nullPtr) {
      for (int i = 0; i < size_; ++i) {
        stat = coeff_[i]*(*xstat)[i];
      }
    }
    return stat;
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real pf(0), one(1);
    Real val = computeValue(obj,x,tol);
    for (int i = 0; i < size_; i++) {
      pf    = plusFunction_->evaluate(val-xstat[i],0);
      val_ += weight_*coeff_[i]/(one-prob_[i])*pf;
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real cvar(0);
    sampler.sumAll(&val_,&cvar,1);
    for (int i = 0; i < size_; i++) {
      cvar += coeff_[i]*xstat[i];
    }
    return cvar;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real pf(0), c(0), one(1);
    Real val = computeValue(obj,x,tol);
    for (int i = 0; i < size_; i++) {
      pf = plusFunction_->evaluate(val-xstat[i],1);
      c  = weight_*coeff_[i]/(one-prob_[i])*pf;
      if (std::abs(c) >= ROL_EPSILON<Real>()) {
        vec_[i] -= c;
        computeGradient(*dualVector_,obj,x,tol);
        g_->axpy(c,*dualVector_);
      }
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    sampler.sumAll(&vec_[0],&gstat[0],size_);
    for (int i = 0; i < size_; i++) {
      gstat[i] += coeff_[i];
    }
    sampler.sumAll(*g_,g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real pf1(0), pf2(0), c(0), one(1);
    Real val = computeValue(obj,x,tol);
    for (int i = 0; i < size_; i++) {
      pf1 = plusFunction_->evaluate(val-xstat[i],1);
      pf2 = plusFunction_->evaluate(val-xstat[i],2);
      if (std::abs(pf2) >= ROL_EPSILON<Real>()) {
        Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
        c        = weight_*coeff_[i]/(one-prob_[i])*pf2*(gv-vstat[i]);
        vec_[i] -= c;
        hv_->axpy(c,*dualVector_);
      }
      if (std::abs(pf1) >= ROL_EPSILON<Real>()) {
        c = weight_*coeff_[i]/(one-prob_[i])*pf1;
        computeHessVec(*dualVector_,obj,v,x,tol);
        RandVarFunctional<Real>::hv_->axpy(c,*dualVector_);
      }
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(&vec_[0],&hvstat[0],size_);
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
