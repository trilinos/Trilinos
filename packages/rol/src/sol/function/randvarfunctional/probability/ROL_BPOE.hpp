// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BPOE_HPP
#define ROL_BPOE_HPP

#include "ROL_RandVarFunctional.hpp"

/** @ingroup stochastic_group 
    \class ROL::BPOE
    \brief Provides the implementation of the buffered probability of exceedance.

    Let \f$(\Omega,\mathcal{F},\mathbb{P})\f$ be a complete space.
    Here, \f$\Omega\f$ is the set of outcomes,
    \f$\mathcal{F}\subseteq 2^\Omega\f$ is a \f$\sigma\f$-algebra of events and
    \f$\mathbb{P}:\mathcal{F}\to[0,1]\f$ is a probability measure.  Moreover,
    let \f$\mathcal{X}\f$ be a class of random variables.

    ROL's BPOE class inherits from ROL::RandVarFunctional which is written in a way
    to exploit parallel sampling.
*/

namespace ROL {

template<class Real>
class BPOE : public RandVarFunctional<Real> {
private:
  Real threshold_;
  Real order_;

  std::vector<Real> hvec_;
  ROL::Ptr<Vector<Real> > dualVec1_, dualVec2_;

  bool firstResetBPOE_;

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

public:
  BPOE(const Real threshold, const Real order=1)
    : RandVarFunctional<Real>(), threshold_(threshold), order_(order), firstResetBPOE_(true) {
    hvec_.resize(5);
  }

  BPOE(ROL::ParameterList &parlist) : RandVarFunctional<Real>(), firstResetBPOE_(true) {
    ROL::ParameterList &list = parlist.sublist("SOL").sublist("Probability").sublist("bPOE");
    threshold_ = list.get<Real>("Threshold");
    order_     = list.get<Real>("Moment Order");
    hvec_.resize(5);
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    if ( firstResetBPOE_ ) {
      dualVec1_ = x.dual().clone();
      dualVec2_ = x.dual().clone();
      firstResetBPOE_ = false;
    }
    dualVec1_->zero();
    dualVec2_->zero();
    hvec_.assign(5,0);
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    const Real zero(0), one(1);
    Real val = computeValue(obj,x,tol);
    Real bp  = xstat[0]*(val-threshold_)+one;
    if ( bp > zero ) {
      val_ += weight_*((order_==one) ? bp : std::pow(bp,order_));
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    const Real one(1);
    Real bpoe(0);
    sampler.sumAll(&val_,&bpoe,1);
    return ((order_==one) ? bpoe : std::pow(bpoe,one/order_));
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real zero(0), one(1), two(2);
    Real val = computeValue(obj,x,tol);
    Real bp  = xstat[0]*(val-threshold_)+one;
    if ( bp > zero ) {
      computeGradient(*dualVector_,obj,x,tol);
      Real pvalp0 = ((order_==one) ? bp : std::pow(bp,order_));
      Real pvalp1 = ((order_==one) ? one : ((order_==two) ? bp : std::pow(bp,order_-one)));
      val_ += weight_ * pvalp0;
      gv_  += weight_ * pvalp1 * (val - threshold_);
      g_->axpy(weight_ * pvalp1, *dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> myvals(2), gvals(2);
    myvals[0] = val_; myvals[1] = gv_;
    sampler.sumAll(&myvals[0],&gvals[0],2);
    if ( gvals[0] > zero) {
      sampler.sumAll(*g_,g);
      Real norm = std::pow(gvals[0],(order_-one)/order_);
      g.scale(xstat[0]/norm);
      gstat[0] = gvals[1]/norm;
    }
    else {
      g.zero();
      gstat[0] = zero;
    }
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real zero(0), one(1), two(2), three(3);
    Real val = computeValue(obj,x,tol);
    Real bp  = xstat[0]*(val-threshold_)+one;
    if ( bp > zero ) {
      // Gradient only
      Real gv     = computeGradVec(*dualVector_,obj,v,x,tol);
      Real pvalp0 = ((order_==one) ? bp : std::pow(bp,order_));
      Real pvalp1 = ((order_==one) ? one
                      : ((order_==two) ? bp : std::pow(bp,order_-one)));
      Real pvalp2 = ((order_==one) ? zero
                      : ((order_==two) ? one
                      : ((order_==three) ? bp : std::pow(bp,order_-two))));
      hvec_[0] += weight_ * pvalp0;
      hvec_[1] += weight_ * pvalp1 * (val-threshold_);
      hvec_[2] += weight_ * pvalp2 * (val-threshold_) * (val-threshold_);
      hvec_[3] += weight_ * pvalp1 * gv; 
      hvec_[4] += weight_ * pvalp2 * (val-threshold_) * gv;
      g_->axpy(weight_ * pvalp1, *dualVector_);
      dualVec1_->axpy(weight_ * pvalp2 * (val-threshold_), *dualVector_);
      dualVec2_->axpy(weight_ * pvalp2 * gv, *dualVector_);
      // Hessian only
      computeHessVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_ * pvalp1, *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1), two(2);
    std::vector<Real> gvals(5);
    sampler.sumAll(&hvec_[0],&gvals[0],5);

    if ( gvals[0] > zero ) {
      Real norm0 = ((order_==one) ? one
                     : ((order_==two) ? std::sqrt(gvals[0])
                       : std::pow(gvals[0],(order_-one)/order_)));
      Real norm1 = ((order_==one) ? gvals[0]
                     : std::pow(gvals[0],(two*order_-one)/order_));
      hvstat[0]  = (order_-one)*((gvals[2]/norm0 - gvals[1]*gvals[1]/norm1)*vstat[0]
                          +xstat[0]*(gvals[4]/norm0 - gvals[3]*gvals[1]/norm1))
                          +(gvals[3]/norm0);

      sampler.sumAll(*hv_,hv);
      hv.scale(xstat[0]/norm0);

      sampler.sumAll(*g_,*hv_);
      Real coeff = -(order_-one)*xstat[0]*(xstat[0]*gvals[3]+vstat[0]*gvals[1])/norm1+vstat[0]/norm0;
      hv.axpy(coeff,*hv_);

      sampler.sumAll(*dualVec1_,*hv_);
      hv.axpy((order_-one)*vstat[0]*xstat[0]/norm0,*hv_);

      sampler.sumAll(*dualVec2_,*hv_);
      hv.axpy((order_-one)*xstat[0]*xstat[0]/norm0,*hv_);
    }
  }
};

}

#endif
