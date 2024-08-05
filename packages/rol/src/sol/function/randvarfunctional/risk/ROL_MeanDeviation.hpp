// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANDEVIATION_HPP
#define ROL_MEANDEVIATION_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

/** @ingroup risk_group
    \class ROL::MeanDeviation
    \brief Provides an interface for the mean plus a sum of arbitrary order
    deviations.

    The mean plus deviations risk measure is
    \f[
       \mathcal{R}(X) = \mathbb{E}[X]
        + \sum_{k=1}^n c_k \mathbb{E}[\wp(X-\mathbb{E}[X])^{p_k}]^{1/p_k}
    \f]
    where \f$\wp:\mathbb{R}\to[0,\infty)\f$ is either the absolute value
    or \f$(x)_+ = \max\{0,x\}\f$, \f$c_k > 0\f$ and \f$p_k\in\mathbb{N}\f$.
    In general, \f$\mathcal{R}\f$ is law-invariant, but not coherent.
    In the specific case that \f$\wp(x) = (x)_+\f$ and \f$c_k\in[0,1]\f$,
    \f$\mathcal{R}\f$ is coherent.  On the other hand,
    the common mean-plus-standard-deviation risk measure (i.e.,
    \f$\wp(x) = |x|\f$, \f$n=1\f$ and \f$p_1 = 2\f$) is not coherent since
    it violates monotonicity.

    When using derivative-based optimization, the user can
    provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PositiveFunction class.
*/

namespace ROL {

template<class Real>
class MeanDeviation : public RandVarFunctional<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Ptr<PositiveFunction<Real> > positiveFunction_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> dev0_;
  std::vector<Real> dev1_;
  std::vector<Real> dev2_;
  std::vector<Real> dev3_;
  std::vector<Real> des0_;
  std::vector<Real> des1_;
  std::vector<Real> des2_;
  std::vector<Real> des3_;
  std::vector<Real> devp_;
  std::vector<Real> gvp1_;
  std::vector<Real> gvp2_;
  std::vector<Real> gvp3_;
  std::vector<Real> gvs1_;
  std::vector<Real> gvs2_;
  std::vector<Real> gvs3_;

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
    NumMoments_ = order_.size();

    dev0_.clear(); dev1_.clear(); dev2_.clear(); dev3_.clear();
    des0_.clear(); des1_.clear(); des2_.clear(); des3_.clear();
    devp_.clear();
    gvp1_.clear(); gvp2_.clear(); gvp3_.clear();
    gvs1_.clear(); gvs2_.clear(); gvs3_.clear();

    dev0_.resize(NumMoments_); dev1_.resize(NumMoments_);
    dev2_.resize(NumMoments_); dev3_.resize(NumMoments_);
    des0_.resize(NumMoments_); des1_.resize(NumMoments_);
    des2_.resize(NumMoments_); des3_.resize(NumMoments_);
    devp_.resize(NumMoments_);
    gvp1_.resize(NumMoments_); gvp2_.resize(NumMoments_);
    gvp3_.resize(NumMoments_);
    gvs1_.resize(NumMoments_); gvs2_.resize(NumMoments_);
    gvs3_.resize(NumMoments_);

    values_    = makePtr<ScalarController<Real>>();
    gradvecs_  = makePtr<ScalarController<Real>>();
    gradients_ = makePtr<VectorController<Real>>();
    hessvecs_  = makePtr<VectorController<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  void clear(void) {
    const Real zero(0);
    dev0_.assign(NumMoments_,zero); dev1_.assign(NumMoments_,zero);
    dev2_.assign(NumMoments_,zero); dev3_.assign(NumMoments_,zero);
    des0_.assign(NumMoments_,zero); des1_.assign(NumMoments_,zero);
    des2_.assign(NumMoments_,zero); des3_.assign(NumMoments_,zero);
    devp_.assign(NumMoments_,zero);
    gvp1_.assign(NumMoments_,zero); gvp2_.assign(NumMoments_,zero);
    gvp3_.assign(NumMoments_,zero);
    gvs1_.assign(NumMoments_,zero); gvs2_.assign(NumMoments_,zero);
    gvs3_.assign(NumMoments_,zero);

    gradvecs_->reset();
    hessvecs_->reset();
  }

  void checkInputs(void) {
    int oSize = order_.size(), cSize = coeff_.size();
    ROL_TEST_FOR_EXCEPTION((oSize!=cSize),std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviation): Order and coefficient arrays have different sizes!");
    Real zero(0), two(2);
    for (int i = 0; i < oSize; i++) {
      ROL_TEST_FOR_EXCEPTION((order_[i] < two), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviation): Element of order array out of range!");
      ROL_TEST_FOR_EXCEPTION((coeff_[i] < zero), std::invalid_argument,
        ">>> ERROR (ROL::MeanDeviation): Element of coefficient array out of range!");
    }
    ROL_TEST_FOR_EXCEPTION(positiveFunction_ == nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::MeanDeviation): PositiveFunction pointer is null!");
    initializeStorage();
  }

public:
  /** \brief Constructor.

      @param[in]     order   is the deviation order
      @param[in]     coeff   is the weight for deviation term
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation risk measure
      with a single deviation.
  */
  MeanDeviation( const Real order, const Real coeff,
                 const Ptr<PositiveFunction<Real> > &pf )
    : RandVarFunctional<Real>(), positiveFunction_(pf) {
    order_.clear(); order_.push_back(order);
    coeff_.clear(); coeff_.push_back(coeff);
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     order   is a vector of deviation orders
      @param[in]     coeff   is a vector of weights for the deviation terms
      @param[in]     pf      is the plus function or an approximation

      This constructor produces a mean plus deviation risk measure
      with an arbitrary number of deviations.
  */
  MeanDeviation( const std::vector<Real> &order,
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
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Mean Plus Deviation" and
      within the "Mean Plus Deviation" sublist should have the following parameters
      \li "Orders" (array of unsigned integers)
      \li "Coefficients" (array of positive scalars)
      \li "Deviation Type" (eighter "Upper" or "Absolute")
      \li A sublist for positive function information.
  */
  MeanDeviation( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation");

    // Get data from parameter list
    order_ = ROL::getArrayFromStringParameter<double>(list,"Orders");
 
    coeff_  = ROL::getArrayFromStringParameter<double>(list,"Coefficients");

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
    val_ += weight_ * val;
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
    //g_->axpy(weight_,*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    //hv_->axpy(weight_,hv);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    // Compute expected value
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    // Compute deviation
    Real diff(0), pf0(0), dev(0), one(1), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0    = positiveFunction_->evaluate(diff,0);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += std::pow(pf0,order_[p]) * weight;
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev += coeff_[p]*std::pow(des0_[p],one/order_[p]);
    }
    // Return mean plus deviation
    return ev + dev;
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
    Real diff(0), pf0(0), pf1(0), c(0), one(1), zero(0), weight(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += weight * std::pow(pf0,order_[p]);
        dev1_[p] += weight * std::pow(pf0,order_[p]-one) * pf1;
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    sampler.sumAll(&dev1_[0],&des1_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev0_[p] = std::pow(des0_[p],one-one/order_[p]);
    }
    // Compute derivative
    dualVector_->zero();
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      c    = zero;
      for ( uint p = 0; p < NumMoments_; p++ ) {
        if ( dev0_[p] > zero ) {
          c += coeff_[p]/dev0_[p] * (std::pow(pf0,order_[p]-one)*pf1 - des1_[p]);
        }
      }
      gradients_->get(*dualVector_,sampler.getMyPoint(i));
      g_->axpy(weight*c,*dualVector_);
    }
    sampler.sumAll(*g_,g);
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
    Real ev = val[0], egv = val[1];
    // Compute deviation
    Real diff(0), pf0(0), pf1(0), pf2(0), zero(0), one(1), two(2);
    Real cg(0), ch(0), diff1(0), diff2(0), diff3(0), weight(0), gv(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += weight * std::pow(pf0,order_[p]);
        dev1_[p] += weight * std::pow(pf0,order_[p]-one) * pf1;
        dev2_[p] += weight * std::pow(pf0,order_[p]-two) * pf1 * pf1;
        dev3_[p] += weight * std::pow(pf0,order_[p]-one) * pf2;
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    sampler.sumAll(&dev1_[0],&des1_[0],NumMoments_);
    sampler.sumAll(&dev2_[0],&des2_[0],NumMoments_);
    sampler.sumAll(&dev3_[0],&des3_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      devp_[p] = std::pow(des0_[p],two-one/order_[p]);
      dev0_[p] = std::pow(des0_[p],one-one/order_[p]);
    }
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      gradvecs_->get(gv,sampler.getMyPoint(i));
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        gvp1_[p] += weight * (std::pow(pf0,order_[p]-one)*pf1-des1_[p]) *
                     (gv - egv);
        gvp2_[p] += weight * (std::pow(pf0,order_[p]-two)*pf1*pf1-des2_[p]) *
                     (gv - egv);
        gvp3_[p] += weight * (std::pow(pf0,order_[p]-one)*pf2-des3_[p]) *
                     (gv - egv);
      }
    }
    sampler.sumAll(&gvp1_[0],&gvs1_[0],NumMoments_);
    sampler.sumAll(&gvp2_[0],&gvs2_[0],NumMoments_);
    sampler.sumAll(&gvp3_[0],&gvs3_[0],NumMoments_);
    // Compute derivative
    dualVector_->zero();
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(diff,sampler.getMyPoint(i));
      weight = sampler.getMyWeight(i);
      diff  -= ev;
      gradvecs_->get(gv,sampler.getMyPoint(i));
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      cg   = one;
      ch   = zero;
      for ( uint p = 0; p < NumMoments_; p++ ) {
        if ( dev0_[p] > zero ) {
          diff1 = std::pow(pf0,order_[p]-one)*pf1-des1_[p];
          diff2 = std::pow(pf0,order_[p]-two)*pf1*pf1*(gv-egv)-gvs2_[p];
          diff3 = std::pow(pf0,order_[p]-one)*pf2*(gv-egv)-gvs3_[p];
          cg   += coeff_[p]*diff1/dev0_[p];
          ch   += coeff_[p]*(((order_[p]-one)*diff2+diff3)/dev0_[p] -
                    (order_[p]-one)*gvs1_[p]*diff1/devp_[p]);
        }
      }
      gradients_->get(*g_,sampler.getMyPoint(i));
      dualVector_->axpy(weight*ch,*g_);
      hessvecs_->get(*hv_,sampler.getMyPoint(i));
      dualVector_->axpy(weight*cg,*hv_);
    }
    sampler.sumAll(*dualVector_,hv);
  }
};

}

#endif
