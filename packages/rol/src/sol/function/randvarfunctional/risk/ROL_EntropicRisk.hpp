// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EXPUTILITY_HPP
#define ROL_EXPUTILITY_HPP

#include "ROL_RandVarFunctional.hpp"

namespace ROL {

/** @ingroup risk_group
    \class ROL::EntropicRisk
    \brief Provides an interface for the entropic risk.

    The entropic risk measure (also called the exponential utility and the
    log-exponential risk measure) is
    \f[
       \mathcal{R}(X) = \lambda
       \log\mathbb{E}\left[\exp\left(\frac{X}{\lambda}\right)\right]
    \f]
    for \f$\lambda > 0\f$.  The entropic risk is convex, translation
    equivariant and monotonic.
*/

template<class Real>
class EntropicRisk : public RandVarFunctional<Real> {
private:
  Real coeff_;

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

  void checkInputs(void) const {
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::EntropicRisk): Rate must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     coeff    is the scale parameter \f$\lambda\f$
  */
  EntropicRisk(const Real coeff = 1)
    : RandVarFunctional<Real>(), coeff_(coeff) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measures"->"Entropic Risk"
      and withing the "Entropic Risk" sublist should have
      \li "Rate" (greater than 0).
  */
  EntropicRisk(ROL::ParameterList &parlist)
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Entropic Risk");
    coeff_ = list.get<Real>("Rate");
    checkInputs();
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * std::exp(coeff_*val);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    return std::log(ev)/coeff_;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = std::exp(coeff_*val);
    val_    += weight_ * ev;
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_*ev,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    Real ev(0), one(1);
    sampler.sumAll(&val_,&ev,1);
    sampler.sumAll(*g_,g);
    g.scale(one/ev);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = std::exp(coeff_*val);
    val_    += weight_ * ev;
    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    gv_     -= weight_ * ev * gv;
    g_->axpy(weight_*ev,*dualVector_);
    hv_->axpy(weight_*coeff_*ev*gv,*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_*ev,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    Real one(1);
    std::vector<Real> myval(2), val(2);
    myval[0] = val_;
    myval[1] = gv_;
    sampler.sumAll(&myval[0],&val[0],2);

    sampler.sumAll(*hv_,hv);
    hv.scale(one/val[0]);

    dualVector_->zero();
    sampler.sumAll(*g_,*dualVector_);
    hv.axpy(coeff_*val[1]/(val[0]*val[0]),*dualVector_);
  }
};

}

#endif
