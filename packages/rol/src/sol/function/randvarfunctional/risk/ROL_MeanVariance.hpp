// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVARIANCE_HPP
#define ROL_MEANVARIANCE_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "ROL_ParameterList.hpp"

/** @ingroup risk_group
    \class ROL::MeanVariance
    \brief Provides an interface for the mean plus a sum of arbitrary order
    variances.

    The mean plus variances risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + \sum_{k=1}^n c_k \mathbb{E}[\wp(X-\mathbb{E}[X])^{p_k}]
    \f]
    where \f$\wp:\mathbb{R}\to[0,\infty)\f$ is either the absolute value
    or \f$(x)_+ = \max\{0,x\}\f$, \f$c_k > 0\f$ and \f$p_k\in\mathbb{N}\f$.
    \f$\mathcal{R}\f$ is law-invariant, but not coherent since it
    violates positive homogeneity.  When \f$\wp(x) = |x|\f$, \f$\mathcal{R}\f$
    also violates monotonicity.

    When using derivative-based optimization, the user can
    provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PositiveFunction class.
*/

namespace ROL {

template<class Real>
class MeanVariance : public RandVarFunctional<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Ptr<PositiveFunction<Real> > positiveFunction_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

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

  void initializeMV(void) {
    values_    = makePtr<ScalarController<Real>>();
    gradvecs_  = makePtr<ScalarController<Real>>();
    gradients_ = makePtr<VectorController<Real>>();
    hessvecs_  = makePtr<VectorController<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  void checkInputs(void) {
    int oSize = order_.size(), cSize = coeff_.size();
    ROL_TEST_FOR_EXCEPTION((oSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanVariance): Order and coefficient arrays have different sizes!");
    Real zero(0), two(2);
    for (int i = 0; i < oSize; i++) {
      ROL_TEST_FOR_EXCEPTION((order_[i] < two), std::invalid_argument,
        ">>> ERROR (ROL::MeanVariance): Element of order array out of range!");
      ROL_TEST_FOR_EXCEPTION((coeff_[i] < zero), std::invalid_argument,
        ">>> ERROR (ROL::MeanVariance): Element of coefficient array out of range!");
    }
    ROL_TEST_FOR_EXCEPTION(positiveFunction_ == nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanVariance): PositiveFunction pointer is null!");
    initializeMV();
  }

public:
  /** \brief Constructor.

      @param[in]     order   is the variance order
      @param[in]     coeff   is the weight for variance term
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus variance risk measure
      with a single variance.
  */
  MeanVariance( const Real order, const Real coeff,
                const Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf) {
    order_.clear(); order_.push_back(order);
    coeff_.clear(); coeff_.push_back(coeff);
    checkInputs();
    NumMoments_ = order_.size();
  }

  /** \brief Constructor.

      @param[in]     order   is a vector of variance orders
      @param[in]     coeff   is a vector of weights for the variance terms
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus variance risk measure
      with an arbitrary number of variances.
  */
  MeanVariance( const std::vector<Real> &order,
                const std::vector<Real> &coeff, 
                const Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf) {
    order_.clear(); coeff_.clear();
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

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Variance" and
      within the "Mean Plus Variance" sublist should have the following parameters
      \li "Orders" (array of unsigned integers)
      \li "Coefficients" (array of positive scalars)
      \li "Deviation Type" (eighter "Upper" or "Absolute")
      \li A sublist for positive function information.
  */
  MeanVariance( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance");
    // Get data from parameter list
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
        ">>> (ROL::MeanVariance): Variance type is not recoginized!");
    }
    // Check inputs
    checkInputs();
    NumMoments_ = order_.size();
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
    Real val(0), diff(0), pf0(0), var(0), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0    = positiveFunction_->evaluate(diff,0);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        val += coeff_[p] * std::pow(pf0,order_[p]) * weight;
      }
    }
    sampler.sumAll(&val,&var,1);
    // Return mean plus deviation
    return ev + var;
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
    Real ev(0), zero(0), one(1);
    sampler.sumAll(&val_,&ev,1);
    sampler.sumAll(*g_,g);
    // Compute deviation
    g_->zero(); dualVector_->zero();
    Real diff(0), pf0(0), pf1(0), c(0), ec(0), ecs(0), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0    = positiveFunction_->evaluate(diff,0);
      pf1    = positiveFunction_->evaluate(diff,1);
      c      = zero;
      for ( uint p = 0; p < NumMoments_; p++ ) {
        c += coeff_[p]*order_[p]*std::pow(pf0,order_[p]-one)*pf1;
      }
      ec += weight*c;
      gradients_->get(*dualVector_,sampler.getMyPoint(i));
      g_->axpy(weight*c,*dualVector_);
    }
    dualVector_->zero();
    sampler.sumAll(&ec,&ecs,1);
    g.scale(one-ecs);
    sampler.sumAll(*g_,*dualVector_);
    g.plus(*dualVector_);
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
    g_->axpy(weight_,*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    // Compute expected value
    std::vector<Real> myval(2), val(2);
    myval[0] = val_;
    myval[1] = gv_;
    sampler.sumAll(&myval[0],&val[0],2);
    Real ev = myval[0], egv = myval[1];
    dualVector_->zero();
    sampler.sumAll(*g_,*dualVector_);
    sampler.sumAll(*hv_,hv);
    // Compute deviation
    g_->zero(); hv_->zero();
    Real diff(0), pf0(0), pf1(0), pf2(0), zero(0), one(1), two(2);
    Real cg(0), ecg(0), ecgs(0), ch(0), ech(0), echs(0), weight(0), gv(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      gradvecs_->get(gv,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0    = positiveFunction_->evaluate(diff,0);
      pf1    = positiveFunction_->evaluate(diff,1);
      pf2    = positiveFunction_->evaluate(diff,2);
      cg     = zero;
      ch     = zero;
      for ( uint p = 0; p < NumMoments_; p++ ) {
        cg += coeff_[p]*order_[p]*(gv-egv)*
                ((order_[p]-one)*std::pow(pf0,order_[p]-two)*pf1*pf1+
                std::pow(pf0,order_[p]-one)*pf2);
        ch += coeff_[p]*order_[p]*std::pow(pf0,order_[p]-one)*pf1;
      }
      ecg += weight*cg;
      ech += weight*ch;
      gradients_->get(*hv_,sampler.getMyPoint(i));
      g_->axpy(weight*cg,*hv_);
      hessvecs_->get(*hv_,sampler.getMyPoint(i));
      g_->axpy(weight*ch,*hv_);
    }
    sampler.sumAll(&ech,&echs,1);
    hv.scale(one-echs);
    sampler.sumAll(&ecg,&ecgs,1);
    hv.axpy(-ecgs,*dualVector_);
    dualVector_->zero();
    sampler.sumAll(*g_,*dualVector_);
    hv.plus(*dualVector_);
  }
};

}

#endif
