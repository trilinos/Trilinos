// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANSEMIDEVIATION_HPP
#define ROL_MEANSEMIDEVIATION_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PlusFunction.hpp"

/** @ingroup risk_group
    \class ROL::MeanSemiDeviation
    \brief Provides an interface for the mean plus upper semideviation of
           order 1.

    The mean plus upper semideviation of order 1 with constant \f$0 < c < 1\f$
    is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
         + c \mathbb{E}\left[(X-\mathbb{E}[X])_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.
    \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class MeanSemiDeviation : public RandVarFunctional<Real> {
private:
  Ptr<PlusFunction<Real> > plusFunction_;
  Real coeff_;

  Ptr<ScalarController<Real>> values_;
  Ptr<ScalarController<Real>> gradvecs_;
  Ptr<VectorController<Real>> gradients_;
  Ptr<VectorController<Real>> hessvecs_;

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
    const Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ < zero), std::invalid_argument,
      ">>> ERROR (ROL::MeanSemiDeviation): Coefficient must be positive!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanSemiDeviation): PlusFunction pointer is null!");
    initializeStorage();
  }

public:

  /** \brief Constructor.

      @param[in]     coeff   is the coefficient scaling the semideviation
      @param[in]     pf      is the plus function or an approximation
  */
  MeanSemiDeviation( const Real coeff, const Ptr<PlusFunction<Real> > &pf )
    : RandVarFunctional<Real>(), plusFunction_(pf), coeff_(coeff) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Semi-Deviation" and
      within the "Mean Plus Semi-Deviation" sublist should have the following parameters
      \li "Coefficient" (between 0 and 1)
      \li A sublist for plus function information.
  */
  MeanSemiDeviation( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Semi-Deviation");
    // Check CVaR inputs
    coeff_ = list.get<Real>("Coefficient");
    // Build (approximate) plus function
    plusFunction_ = makePtr<PlusFunction<Real>>(list);
    // Check Inputs
    checkInputs();
  }

  void setStorage(const Ptr<ScalarController<Real>> &value_storage,
                  const Ptr<VectorController<Real>> &gradient_storage) {
    values_    = value_storage;
    gradients_ = gradient_storage;
    RandVarFunctional<Real>::setStorage(values_,gradients_);
  }

  void setHessVecStorage(const Ptr<ScalarController<Real>> &gradvec_storage,
                         const Ptr<VectorController<Real>> &hessvec_storage) {
    gradvecs_ = gradvec_storage;
    hessvecs_ = hessvec_storage;
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }
 
  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    clear();
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * val;
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    // Compute expected value
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    // Compute deviation
    Real diff(0), pf(0), dev(0), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf    += weight * plusFunction_->evaluate(diff,0);
    }
    sampler.sumAll(&pf,&dev,1);
    // Return mean plus deviation
    return ev + coeff_ * dev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * val;
    computeGradient(*dualVector_,obj,x,tol);
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
    Real diff(0), dev(0), pf (0), c(0), one(1), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf    += weight * plusFunction_->evaluate(diff,1);
    }
    sampler.sumAll(&pf,&dev,1);
    // Compute derivative
    g_->zero(); dualVector_->zero();
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf     = plusFunction_->evaluate(diff,1);
      c      = one + coeff_ * (pf - dev);
      gradients_->get(*dualVector_, sampler.getMyPoint(i));
      g_->axpy(weight * c, *dualVector_);
    }
    sampler.sumAll(*g_, g);
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
    computeHessVec(*dualVector_,obj,v,x,tol);
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
    std::vector<Real> mval = {val_, gv_};
    std::vector<Real> gval(2,0);
    sampler.sumAll(&mval[0],&gval[0],2);
    Real ev = gval[0], egv = gval[1];
    // Compute deviation
    Real diff(0), pf1(0), pf2(0), weight(0), gv(0), c(0);
    std::vector<Real> pf(2,0), dev(2,0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff, sampler.getMyPoint(i));
      gradvecs_->get(gv, sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf[0] += weight * plusFunction_->evaluate(diff,1);
      pf[1] += weight * plusFunction_->evaluate(diff,2) * (gv - egv);
    }
    sampler.sumAll(&pf[0],&dev[0],2);
    hv_->zero(); dualVector_->zero();
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff, sampler.getMyPoint(i));
      gradvecs_->get(gv, sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf1    = plusFunction_->evaluate(diff,1);
      c      = one + coeff_ * (pf1 - dev[0]);
      hessvecs_->get(*dualVector_, sampler.getMyPoint(i));
      hv_->axpy(weight * c, *dualVector_);
      pf2    = plusFunction_->evaluate(diff,2) * (gv - egv);
      c      = coeff_ * (pf2 - dev[1]);
      gradients_->get(*dualVector_, sampler.getMyPoint(i));
      hv_->axpy(weight * c, *dualVector_);
    }
    sampler.sumAll(*hv_, hv);
  }
};

}

#endif
