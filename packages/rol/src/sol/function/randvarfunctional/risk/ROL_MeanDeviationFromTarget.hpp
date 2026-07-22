// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANDEVIATIONFROMTARGET_HPP
#define ROL_MEANDEVIATIONFROMTARGET_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

/** @ingroup risk_group
    \class ROL::MeanDeviationFromTarget
    \brief Provides an interface for the mean plus a sum of arbitrary order
    deviations from targets.

    The mean plus deviations from targets risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + \sum_{k=1}^n c_k \mathbb{E}[\wp(X-t_k)^{p_k}]^{1/p_k}
    \f]
    where \f$\wp:\mathbb{R}\to[0,\infty)\f$ is either the absolute value
    or \f$(x)_+ = \max\{0,x\}\f$, \f$c_k > 0\f$ and \f$p_k\in\mathbb{N}\f$.
    In general, \f$\mathcal{R}\f$ is law-invariant, but not coherent since
    it violates translation equivariance.

    When using derivative-based optimization, the user can
    provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PositiveFunction class.
*/

namespace ROL {

template<class Real>
class MeanDeviationFromTarget : public RandVarFunctional<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Ptr<PositiveFunction<Real> > positiveFunction_;
  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> pval_; 
  std::vector<Real> pgv_; 

  std::vector<Ptr<Vector<Real> > > pg0_;
  std::vector<Ptr<Vector<Real> > > pg_;
  std::vector<Ptr<Vector<Real> > > phv_;

  bool firstResetMDT_;

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

  void initializeMDT(void) {
    // Initialize additional storage
    pg_.clear(); pg0_.clear(); phv_.clear(); pval_.clear(); pgv_.clear();
    pg_.resize(NumMoments_);
    pg0_.resize(NumMoments_);
    phv_.resize(NumMoments_);
    pval_.resize(NumMoments_);
    pgv_.resize(NumMoments_);
  }

  void checkInputs(void) {
    int oSize = order_.size(), cSize = coeff_.size(), tSize = target_.size();
    ROL_TEST_FOR_EXCEPTION((oSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): Order and coefficient arrays have different sizes!");
    ROL_TEST_FOR_EXCEPTION((oSize!=tSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): Order and target arrays have different sizes!");
    Real zero(0), two(2);
    for (int i = 0; i < oSize; i++) {
      ROL_TEST_FOR_EXCEPTION((order_[i] < two), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviationFromTarget): Element of order array out of range!");
      ROL_TEST_FOR_EXCEPTION((coeff_[i] < zero), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviationFromTarget): Element of coefficient array out of range!");
    }
    ROL_TEST_FOR_EXCEPTION(positiveFunction_ == nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviationFromTarget): PositiveFunction pointer is null!");
    initializeMDT();
  }

public:
  /** \brief Constructor.

      @param[in]     target  is the target scalar value
      @param[in]     order   is the deviation order
      @param[in]     coeff   is the weight for deviation term
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation from target risk measure
      with a single deviation.
  */
  MeanDeviationFromTarget( const Real target, const Real order, const Real coeff,
                           const Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf), firstResetMDT_(true) {
    order_.clear();  order_.push_back(order);
    coeff_.clear();  coeff_.push_back(coeff);
    target_.clear(); target_.push_back(target);
    NumMoments_ = order_.size();
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     target  is a vector of targets
      @param[in]     order   is a vector of deviation orders
      @param[in]     coeff   is a vector of weights for the deviation terms
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation from target risk measure
      with an arbitrary number of deviations.
  */
  MeanDeviationFromTarget( const std::vector<Real> &target,
                           const std::vector<Real> &order,
                           const std::vector<Real> &coeff, 
                           const Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf), firstResetMDT_(true) {
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
    NumMoments_ = order_.size();
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Deviation From Target" and
      within the "Mean Plus Deviation From Target" sublist should have the following parameters
      \li "Targets" (array of scalars)
      \li "Orders" (array of unsigned integers)
      \li "Coefficients" (array of positive scalars)
      \li "Deviation Type" (eighter "Upper" or "Absolute")
      \li A sublist for positive function information.
  */
  MeanDeviationFromTarget( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>(), firstResetMDT_(true) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation From Target");

    // Get data from parameter list
    target_ = ROL::getArrayFromStringParameter<double>(list,"Targets");

    order_ = ROL::getArrayFromStringParameter<double>(list,"Orders");

    coeff_ = ROL::getArrayFromStringParameter<double>(list,"Coefficients");

    // Build (approximate) positive function
    std::string type = list.get<std::string>("Deviation Type");
    if ( type == "Upper" ) {
      positiveFunction_ = makePtr<PlusFunction<Real>>(list);
    }
    else if ( type == "Absolute" ) {
      positiveFunction_ = makePtr<AbsoluteValue<Real>>(list);
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::MeanDeviation): Deviation type is not recoginized!");
    }
    // Check inputs
    NumMoments_ = order_.size();
    checkInputs();
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    if (firstResetMDT_) {
      for ( uint p = 0; p < NumMoments_; p++ ) {
        pg0_[p] = x.dual().clone();
        pg_[p]  = x.dual().clone();
        phv_[p] = x.dual().clone();
      }
      firstResetMDT_  = false;
    }
    Real zero(0);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      pg0_[p]->zero(); pg_[p]->zero(); phv_[p]->zero();
      pval_[p] = zero; pgv_[p] = zero;
    }
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
      pval_[p] += weight_ * std::pow(pf0,order_[p]);
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    const Real one(1);
    Real dev(0);
    sampler.sumAll(&val_,&dev,1);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&pval_[0],&pval_sum[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev += coeff_[p] * std::pow(pval_sum[p],one/order_[p]);
    }
    return dev;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real diff(0), pf0(0), pf1(0), c(0), one(1);
    Real val = computeValue(obj,x,tol);
    computeGradient(*dualVector_,obj,x,tol);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      c    = std::pow(pf0,order_[p]-one) * pf1;
      (pg_[p])->axpy(weight_ * c,*dualVector_);
      pval_[p] += weight_ * std::pow(pf0,order_[p]);
    }
    g_->axpy(weight_,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1);
    sampler.sumAll(*g_,g);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&pval_[0],&pval_sum[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        dualVector_->zero();
        sampler.sumAll(*(pg_[p]),*dualVector_);
        g.axpy(coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]),*dualVector_);
      }
    }
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real one(1), two(2);
    Real diff(0), pf0(0), pf1(0), pf2(0), p0(0), p1(0), p2(0), c(0);
    Real val = computeValue(obj,x,tol);
    Real gv  = computeGradVec(*g_,obj,v,x,tol);
    computeHessVec(*dualVector_,obj,v,x,tol);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val - target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      pf2 = positiveFunction_->evaluate(diff,2);
      p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-one);
      p2   = std::pow(pf0,order_[p]-two);
      c    = -(order_[p]-one)*p1*pf1;
      pg0_[p]->axpy(weight_*c,*g_);
      c    = gv*((order_[p]-one)*p2*pf1*pf1 + p1*pf2);
      pg_[p]->axpy(weight_*c,*g_);
      c    = p1*pf1;
      phv_[p]->axpy(weight_*c,*dualVector_);
      pval_[p] += weight_*p0;
      pgv_[p]  += weight_*p1*pf1*gv;
    }
    hv_->axpy(weight_,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1), two(2);
    sampler.sumAll(*hv_,hv);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    std::vector<Real> pgv_sum(NumMoments_);
    sampler.sumAll(&(pgv_)[0],&pgv_sum[0],NumMoments_);
    Real c(0);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        dualVector_->zero();
        sampler.sumAll(*(pg0_[p]),*dualVector_);
        c = coeff_[p]*(pgv_sum[p]/std::pow(pval_sum[p],two-one/order_[p]));
        hv.axpy(c,*dualVector_);

        dualVector_->zero();
        sampler.sumAll(*(pg_[p]),*dualVector_);
        c = coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]);
        hv.axpy(c,*dualVector_);

        dualVector_->zero();
        sampler.sumAll(*(phv_[p]),*dualVector_);
        hv.axpy(c,*dualVector_);
      }
    }
  }
};

}

#endif
