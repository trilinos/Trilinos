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

#ifndef ROL_MEANDEVIATIONFROMTARGET_HPP
#define ROL_MEANDEVIATIONFROMTARGET_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace ROL {

template<class Real>
class MeanDeviationFromTarget : public RiskMeasure<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;
  Teuchos::RCP<Vector<Real> > dualVector3_;
  Teuchos::RCP<Vector<Real> > dualVector4_;

  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> pval_; 
  std::vector<Real> pgv_; 

  std::vector<Teuchos::RCP<Vector<Real> > > pg0_;
  std::vector<Teuchos::RCP<Vector<Real> > > pg_;
  std::vector<Teuchos::RCP<Vector<Real> > > phv_;

  bool firstReset_;

  void initialize(void) {
    // Initialize additional storage
    pg_.clear(); pg0_.clear(); phv_.clear(); pval_.clear(); pgv_.clear();
    pg_.resize(NumMoments_);
    pg0_.resize(NumMoments_);
    phv_.resize(NumMoments_);
    pval_.resize(NumMoments_);
    pgv_.resize(NumMoments_);
  }

public:

  MeanDeviationFromTarget( Real target, Real order, Real coeff,
                           Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    Real zero(0), one(1), two(2);
    // Initialize storage for problem data
    target_.clear(); order_.clear(); coeff_.clear();
    target_.push_back(target);
    order_.push_back((order < two) ? two : order);
    coeff_.push_back((coeff < zero) ? one : coeff);
    NumMoments_ = order_.size();
    initialize();
  }

  MeanDeviationFromTarget( std::vector<Real> &target, std::vector<Real> &order, std::vector<Real> &coeff, 
                           Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    Real zero(0), one(1), two(2);
    // Initialize storage for problem data
    NumMoments_ = order.size();
    target_.clear(); order_.clear(); coeff_.clear();
    if ( NumMoments_ != target.size() ) {
      target.resize(NumMoments_,zero);
    }
    if ( NumMoments_ != coeff.size() ) {
      coeff.resize(NumMoments_,one);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      target_.push_back(target[i]);
      order_.push_back((order[i] < two) ? two : order[i]);
      coeff_.push_back((coeff[i] < one) ? one : coeff[i]);
    }
    initialize();
  }

  MeanDeviationFromTarget( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Real zero(0), one(1), two(2);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation From Target");
    // Get data from parameter list
    Teuchos::Array<Real> target
      = Teuchos::getArrayFromStringParameter<double>(list,"Targets");
    Teuchos::Array<Real> order
      = Teuchos::getArrayFromStringParameter<double>(list,"Orders");
    Teuchos::Array<Real> coeff
      = Teuchos::getArrayFromStringParameter<double>(list,"Coefficients");
    // Check inputs
    NumMoments_ = order.size();
    target_.clear(); order_.clear(); coeff_.clear();
    if ( NumMoments_ != static_cast<uint>(target.size()) ) {
      target.resize(NumMoments_,zero);
    }
    if ( NumMoments_ != static_cast<uint>(coeff.size()) ) {
      coeff.resize(NumMoments_,one);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      target_.push_back(target[i]);
      order_.push_back((order[i] < two) ? two : order[i]);
      coeff_.push_back((coeff[i] < zero) ? one : coeff[i]);
    }
    initialize();
    // Build (approximate) positive function
    if ( list.get("Deviation Type","Upper") == "Upper" ) {
      positiveFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    }
    else {
      positiveFunction_ = Teuchos::rcp(new AbsoluteValue<Real>(list));
    }
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    if (firstReset_) {
      for ( uint p = 0; p < NumMoments_; p++ ) {
        pg0_[p] = (x0->dual()).clone();
        pg_[p]  = (x0->dual()).clone();
        phv_[p] = (x0->dual()).clone();
      }
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      dualVector3_ = (x0->dual()).clone();
      dualVector4_ = (x0->dual()).clone();
      firstReset_  = false;
    }
    for ( uint p = 0; p < NumMoments_; p++ ) {
      pg0_[p]->zero(); pg_[p]->zero(); phv_[p]->zero();
      pval_[p] = zero; pgv_[p] = zero;
    }
    dualVector1_->zero(); dualVector2_->zero();
    dualVector3_->zero(); dualVector4_->zero();
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
  }
  
  void update(const Real val, const Real weight) {
    Real diff(0), pf0(0);
    RiskMeasure<Real>::val_ += weight * val;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff(0), pf0(0), pf1(0), c(0), one(1);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      c    = std::pow(pf0,order_[p]-one) * pf1;
      (pg_[p])->axpy(weight * c,g);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
    RiskMeasure<Real>::g_->axpy(weight,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff(0), pf0(0), pf1(0), pf2(0), p0(0), p1(0), p2(0), c(0), one(1), two(2);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val - target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      pf2 = positiveFunction_->evaluate(diff,2);
      p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-one);
      p2   = std::pow(pf0,order_[p]-two);
      c    = -(order_[p]-one)*p1*pf1;
      pg0_[p]->axpy(weight*c,g);
      c    = gv*((order_[p]-one)*p2*pf1*pf1 + p1*pf2);
      pg_[p]->axpy(weight*c,g);
      c    = p1*pf1;
      phv_[p]->axpy(weight*c,hv);
      pval_[p] += weight*p0;
      pgv_[p]  += weight*p1*pf1*gv;
    }
    RiskMeasure<Real>::hv_->axpy(weight,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, dev(0), one(1);
    sampler.sumAll(&val,&dev,1);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev += coeff_[p] * std::pow(pval_sum[p],one/order_[p]);
    }
    return dev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    Teuchos::RCP<Vector<Real> > pg;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        pg = (pg_[p])->clone();
        sampler.sumAll(*(pg_[p]),*pg);
        dualVector1_->axpy(coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]),*pg);
      }
    }
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector1_);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1), two(2);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);
    std::vector<Real> pval_sum(NumMoments_);
    sampler.sumAll(&(pval_)[0],&pval_sum[0],NumMoments_);
    std::vector<Real> pgv_sum(NumMoments_);
    sampler.sumAll(&(pgv_)[0],&pgv_sum[0],NumMoments_);
    Real c(0);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      if ( pval_sum[p] > zero ) {
        sampler.sumAll(*(pg_[p]),*dualVector2_);
        sampler.sumAll(*(pg0_[p]),*dualVector3_);
        sampler.sumAll(*(phv_[p]),*dualVector4_);
        c = coeff_[p]*(pgv_sum[p]/std::pow(pval_sum[p],two-one/order_[p]));
        dualVector1_->axpy(c,*dualVector3_);
        c = coeff_[p]/std::pow(pval_sum[p],one-one/order_[p]);
        dualVector1_->axpy(c,*dualVector2_);
        dualVector1_->axpy(c,*dualVector4_);
      }
    }
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector1_);
  }
};

}

#endif
