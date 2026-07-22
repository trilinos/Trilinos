// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVARIANCEFROMTARGET_HPP
#define ROL_MEANVARIANCEFROMTARGET_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "ROL_ParameterList.hpp"

/** @ingroup risk_group
    \class ROL::MeanVarianceFromTarget
    \brief Provides an interface for the mean plus a sum of arbitrary order
    variances from targets.

    The mean plus variances from targets risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + \sum_{k=1}^n c_k \mathbb{E}[\wp(X-t_k)^{p_k}]
    \f]
    where \f$\wp:\mathbb{R}\to[0,\infty)\f$ is either the absolute value
    or \f$(x)_+ = \max\{0,x\}\f$, \f$c_k > 0\f$ and \f$p_k\in\mathbb{N}\f$.
    \f$\mathcal{R}\f$ is law-invariant, but not coherent since it
    violates positive homogeneity and translation equivariance.

    When using derivative-based optimization, the user can
    provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PositiveFunction class.
*/

namespace ROL {

template<class Real>
class MeanVarianceFromTarget : public RandVarFunctional<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  ROL::Ptr<PositiveFunction<Real> > positiveFunction_;
  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

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
    int oSize = order_.size(), cSize = coeff_.size();
    ROL_TEST_FOR_EXCEPTION((oSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanVarianceFromTarget): Order and coefficient arrays have different sizes!");
    Real zero(0), two(2);
    for (int i = 0; i < oSize; i++) {
      ROL_TEST_FOR_EXCEPTION((order_[i] < two), std::invalid_argument,
        ">>> ERROR (ROL::MeanVarianceFromTarget): Element of order array out of range!");
      ROL_TEST_FOR_EXCEPTION((coeff_[i] < zero), std::invalid_argument,
        ">>> ERROR (ROL::MeanVarianceFromTarget): Element of coefficient array out of range!");
    }
    ROL_TEST_FOR_EXCEPTION(positiveFunction_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanVarianceFromTarget): PositiveFunction pointer is null!");
  }

public:
  /** \brief Constructor.

      @param[in]     target  is the scalar target
      @param[in]     order   is the variance order
      @param[in]     coeff   is the weight for variance term
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus variance from target risk measure
      with a single variance.
  */
  MeanVarianceFromTarget( const Real target, const Real order, const Real coeff,
                          const ROL::Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf) {
    target_.clear(); target_.push_back(target);
    order_.clear();  order_.push_back(order);
    coeff_.clear();  coeff_.push_back(coeff);
    checkInputs();
    NumMoments_ = order_.size();
  }

  /** \brief Constructor.

      @param[in]     target  is a vector of targets
      @param[in]     order   is a vector of variance orders
      @param[in]     coeff   is a vector of weights for the variance terms
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus variance from target risk measure
      with an arbitrary number of variances.
  */
  MeanVarianceFromTarget( const std::vector<Real> &target,
                          const std::vector<Real> &order,
                          const std::vector<Real> &coeff, 
                          const ROL::Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf) {
    target_.clear(); order_.clear(); coeff_.clear();
    for ( uint i = 0; i < target.size(); i++ ) {
      target_.push_back(target[i]);
    }
    for ( uint i = 0; i < order.size(); i++ ) {
      order_.push_back(order[i]);
    }
    for ( uint i = 0; i < coeff.size(); i++ ) {
      coeff_.push_back(coeff[i]);
    }
    checkInputs();
    NumMoments_ = order_.size();
  }
  
  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Variance From Target" and
      within the "Mean Plus Variance From Target" sublist should have the following parameters
      \li "Targets" (array of scalars)
      \li "Orders" (array of unsigned integers)
      \li "Coefficients" (array of positive scalars)
      \li "Deviation Type" (eighter "Upper" or "Absolute")
      \li A sublist for positive function information.
  */
  MeanVarianceFromTarget( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance From Target");
    // Get data from parameter list
    target_ = ROL::getArrayFromStringParameter<double>(list,"Targets");
    order_  = ROL::getArrayFromStringParameter<double>(list,"Orders");
    coeff_  = ROL::getArrayFromStringParameter<double>(list,"Coefficients");

    // Build (approximate) positive function
    std::string type = list.get<std::string>("Deviation Type");
    if ( type == "Upper" ) {
      positiveFunction_ = ROL::makePtr<PlusFunction<Real>>(list);
    }
    else if ( type == "Absolute" ) {
      positiveFunction_ = ROL::makePtr<AbsoluteValue<Real>>(list);
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::MeanDeviation): Deviation type is not recoginized!");
    }
    // Check inputs
    checkInputs();
    NumMoments_ = order_.size();
  }
  
  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real diff(0), pf0(0);
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * val;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      val_ += weight_ * coeff_[p] * std::pow(pf0,order_[p]);
    }
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real diff(0), pf0(0), pf1(0), c(1), one(1);
    Real val = computeValue(obj,x,tol);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      c   += order_[p]*coeff_[p]*std::pow(pf0,order_[p]-one)*pf1;
    }
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_ * c,*dualVector_);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real diff(0), pf0(0), pf1(0), pf2(0), p1(0), p2(0), ch(1), cg(0), one(1), two(2);
    Real val = computeValue(obj,x,tol);
    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val - target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      //p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-one);
      p2   = std::pow(pf0,order_[p]-two);
      cg  += order_[p]*coeff_[p]*gv*( (order_[p]-one)*p2*pf1*pf1 + p1*pf2 );
      ch  += order_[p]*coeff_[p]*p1*pf1;
    }
    hv_->axpy(weight_*cg,*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_*ch,*dualVector_);
  }
};

}

#endif
