// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_BPOE_HPP
#define ROL_BPOE_HPP

#include "ROL_RiskMeasure.hpp"

/** @ingroup stochastic_group 
    \class ROL::BPOE
    \brief Provides the implementation of the buffered probability of exceedance.

    Let \f$(\Omega,\mathcal{F},\mathbb{P})\f$ be a complete space.
    Here, \f$\Omega\f$ is the set of outcomes,
    \f$\mathcal{F}\subseteq 2^\Omega\f$ is a \f$\sigma\f$-algebra of events and
    \f$\mathbb{P}:\mathcal{F}\to[0,1]\f$ is a probability measure.  Moreover,
    let \f$\mathcal{X}\f$ be a class of random variables.

    ROL's BPOE class inherits from ROL::RiskMeasure which is written in a way
    to exploit parallel sampling.
*/

namespace ROL {

template<class Real>
class BPOE : public RiskMeasure<Real> {
private:
  Real threshold_;
  Real order_;

  Real xvar_, vvar_;

  std::vector<Real> hvec_;
  ROL::SharedPointer<Vector<Real> > dualVec1_, dualVec2_;

  bool firstReset_;

public:
  BPOE(const Real threshold, const Real order=1)
    : RiskMeasure<Real>(), threshold_(threshold), order_(order), firstReset_(true) {
    hvec_.resize(5);
  }

  BPOE(ROL::ParameterList &parlist) : RiskMeasure<Real>(), firstReset_(true) {
    threshold_ = parlist.sublist("SOL").sublist("Risk Measure").sublist("bPOE").get("Threshold",1.0);
    order_     = parlist.sublist("SOL").sublist("Risk Measure").sublist("bPOE").get("Moment Order",1.0);
    hvec_.resize(5);
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    xvar_ = (*dynamic_cast<const RiskVector<Real>&>(x).getStatistic(comp,index))[0];
  }

  void reset(ROL::SharedPointer<Vector<Real> > &x0, const Vector<Real> &x,
                     ROL::SharedPointer<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = ROL::constPointerCast<Vector<Real> >(
         dynamic_cast<const RiskVector<Real>&>(v).getVector());
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    vvar_ = (*dynamic_cast<const RiskVector<Real>&>(v).getStatistic(comp,index))[0];

    if ( firstReset_ ) {
      dualVec1_ = (x0->dual()).clone();
      dualVec2_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVec1_->zero();
    dualVec2_->zero();
    hvec_.assign(5,0);
  }

  void update(const Real val, const Real weight) {
    const Real zero(0), one(1);
    Real bp = xvar_*(val-threshold_)+one;
    if ( bp > zero ) {
      RiskMeasure<Real>::val_ += weight
        * ((order_==one) ? bp : std::pow(bp,order_));
    }
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    const Real one(1);
    Real val = RiskMeasure<Real>::val_, bpoe(0);
    sampler.sumAll(&val,&bpoe,1);
    return ((order_==one) ? bpoe : std::pow(bpoe,one/order_));
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    const Real zero(0), one(1), two(2);
    Real bp = xvar_*(val-threshold_)+one;
    if ( bp > zero ) {
      Real pvalp0 = ((order_==one) ? bp : std::pow(bp,order_));
      Real pvalp1 = ((order_==one) ? one : ((order_==two) ? bp : std::pow(bp,order_-one)));
      RiskMeasure<Real>::val_ += weight * pvalp0;
      RiskMeasure<Real>::gv_  += weight * pvalp1 * (val - threshold_);
      RiskMeasure<Real>::g_->axpy(weight * pvalp1, g);
    }
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> myvals(2), gvals(2);
    myvals[0] = RiskMeasure<Real>::val_;
    myvals[1] = RiskMeasure<Real>::gv_;
    sampler.sumAll(&myvals[0],&gvals[0],2);

    Real gvar(0);
    if ( gvals[0] > zero) {
      ROL::SharedPointer<Vector<Real> > gvec
        = dynamic_cast<RiskVector<Real>&>(g).getVector();
      sampler.sumAll(*(RiskMeasure<Real>::g_),*gvec);
      Real norm = std::pow(gvals[0],(order_-one)/order_);
      gvec->scale(xvar_/norm);
      gvar = gvals[1]/norm;
    }
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    dynamic_cast<RiskVector<Real>&>(g).setStatistic(gvar,comp,index);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    const Real zero(0), one(1), two(2), three(3);
    Real bp = xvar_*(val-threshold_)+one;
    if ( bp > zero ) {
      Real pvalp0 = ((order_==one) ? bp : std::pow(bp,order_));
      Real pvalp1 = ((order_==one) ? one
                      : ((order_==two) ? bp : std::pow(bp,order_-one)));
      Real pvalp2 = ((order_==one) ? zero
                      : ((order_==two) ? one
                      : ((order_==three) ? bp : std::pow(bp,order_-two))));
      hvec_[0] += weight * pvalp0;
      hvec_[1] += weight * pvalp1 * (val-threshold_);
      hvec_[2] += weight * pvalp2 * (val-threshold_) * (val-threshold_);
      hvec_[3] += weight * pvalp1 * gv; 
      hvec_[4] += weight * pvalp2 * (val-threshold_) * gv;
      RiskMeasure<Real>::g_->axpy(weight * pvalp1, g);
      dualVec1_->axpy(weight * pvalp2 * (val-threshold_), g);
      dualVec2_->axpy(weight * pvalp2 * gv, g);
      RiskMeasure<Real>::hv_->axpy(weight * pvalp1, hv);
    }
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    const Real zero(0), one(1), two(2);
    std::vector<Real> gvals(5);
    sampler.sumAll(&hvec_[0],&gvals[0],5);

    Real hvar(0);
    if ( gvals[0] > zero ) {
      ROL::SharedPointer<Vector<Real> > hvec
        = dynamic_cast<RiskVector<Real>&>(hv).getVector();
      Real norm0 = ((order_==one) ? one
                     : ((order_==two) ? std::sqrt(gvals[0])
                       : std::pow(gvals[0],(order_-one)/order_)));
      Real norm1 = ((order_==one) ? gvals[0]
                     : std::pow(gvals[0],(two*order_-one)/order_));
      hvar = (order_-one)*((gvals[2]/norm0 - gvals[1]*gvals[1]/norm1)*vvar_
                          +xvar_*(gvals[4]/norm0 - gvals[3]*gvals[1]/norm1))
                          +(gvals[3]/norm0);

      sampler.sumAll(*RiskMeasure<Real>::hv_,*hvec);
      hvec->scale(xvar_/norm0);

      sampler.sumAll(*RiskMeasure<Real>::g_,*RiskMeasure<Real>::hv_);
      Real coeff = -(order_-one)*xvar_*(xvar_*gvals[3]+vvar_*gvals[1])/norm1+vvar_/norm0;
      hvec->axpy(coeff,*RiskMeasure<Real>::hv_);

      sampler.sumAll(*dualVec1_,*RiskMeasure<Real>::hv_);
      hvec->axpy((order_-one)*vvar_*xvar_/norm0,*RiskMeasure<Real>::hv_);

      sampler.sumAll(*dualVec2_,*RiskMeasure<Real>::hv_);
      hvec->axpy((order_-one)*xvar_*xvar_/norm0,*RiskMeasure<Real>::hv_);
    }
    int index = RiskMeasure<Real>::getIndex();
    int comp  = RiskMeasure<Real>::getComponent();
    dynamic_cast<RiskVector<Real>&>(hv).setStatistic(hvar,comp,index);
  }
};

}

#endif
