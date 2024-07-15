// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANSEMIDEVIATIONFROMTARGET_HPP
#define ROL_MEANSEMIDEVIATIONFROMTARGET_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PlusFunction.hpp"

/** @ingroup risk_group
    \class ROL::MeanSemiDeviation
    \brief Provides an interface for the mean plus upper semideviation
           from target of order 1.

    The mean plus upper semideviation from target of order 1 with constant
    \f$0 < c < 1\f$ and target \f$t\in\mathbb{R}\f$ is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
         + c \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.
    \f$\mathcal{R}\f$ is a law-invariant risk measure.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class MeanSemiDeviationFromTarget : public RandVarFunctional<Real> {
private:
  Ptr<PlusFunction<Real> > plusFunction_;
  Real coeff_, target_;

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
    const Real zero(0);
    ROL_TEST_FOR_EXCEPTION((coeff_ < zero), std::invalid_argument,
      ">>> ERROR (ROL::MeanPlusSemiDeviationFromTarget): Coefficient must be positive!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanSemiDeviation): PlusFunction pointer is null!");
  }

public:

  /** \brief Constructor.

      @param[in]     coeff   is the coefficient scaling the semideviation
      @param[in]     target  is the target scalar
      @param[in]     pf      is the plus function or an approximation
  */
  MeanSemiDeviationFromTarget( const Real coeff, const Real target,
                               const Ptr<PlusFunction<Real>> &pf )
    : RandVarFunctional<Real>(), plusFunction_(pf), coeff_(coeff), target_(target) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Semi-Deviation" and
      within the "Mean Plus Semi-Deviation" sublist should have the following parameters
      \li "Coefficient" (between 0 and 1)
      \li "Target"
      \li A sublist for plus function information.
  */
  MeanSemiDeviationFromTarget( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Semi-Deviation From Target");
    // Check inputs
    coeff_  = list.get<Real>("Coefficient");
    target_ = list.get<Real>("Target");
    // Build (approximate) plus function
    plusFunction_ = makePtr<PlusFunction<Real>>(list);
    // Check Inputs
    checkInputs();
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real pf  = plusFunction_->evaluate(val-target_,0);
    val_    += weight_ * (val + coeff_ * pf);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    return ev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real one(1);
    Real val = computeValue(obj,x,tol);
    Real pf  = plusFunction_->evaluate(val-target_,1);
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_ * (one + coeff_ * pf), *dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*g_, g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real one(1);
    Real val = computeValue(obj,x,tol);
    Real pf1 = plusFunction_->evaluate(val-target_,1);
    Real pf2 = plusFunction_->evaluate(val-target_,2);
    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_ * coeff_ * pf2 * gv, *dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_ * (one + coeff_ * pf1), *dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*hv_, hv);
  }
};

}

#endif
