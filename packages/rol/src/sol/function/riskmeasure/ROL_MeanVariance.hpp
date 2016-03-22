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

#ifndef ROL_MEANVARIANCE_HPP
#define ROL_MEANVARIANCE_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PositiveFunction.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace ROL {

template<class Real>
class MeanVariance : public RiskMeasure<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;
  Teuchos::RCP<Vector<Real> > dualVector3_;
  Teuchos::RCP<Vector<Real> > dualVector4_;

  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> weights_;
  std::vector<Real> value_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > gradient_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > hessvec_storage_;
  std::vector<Real> gradvec_storage_;

  bool firstReset_;

public:

  MeanVariance( Real order, Real coeff,
                Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    Real zero(0), one(1), two(2);
    order_.clear(); coeff_.clear();
    order_.push_back((order < two) ? two : order);
    coeff_.push_back((coeff < zero) ? one : coeff);
    NumMoments_ = order_.size();
  }
  MeanVariance( std::vector<Real> &order, std::vector<Real> &coeff, 
                Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    Real zero(0), one(1), two(2);
    order_.clear(); coeff_.clear();
    NumMoments_ = order.size();
    if ( NumMoments_ != coeff.size() ) {
      coeff.resize(NumMoments_,one);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      order_.push_back((order[i] < two) ? two : order[i]);
      coeff_.push_back((coeff[i] < zero) ? one : coeff[i]);
    }
  }
  MeanVariance( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Real zero(0), one(1), two(2);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance");
    // Get data from parameter list
    Teuchos::Array<Real> order
      = Teuchos::getArrayFromStringParameter<double>(list,"Orders");
    Teuchos::Array<Real> coeff
      = Teuchos::getArrayFromStringParameter<double>(list,"Coefficients");
    // Check inputs
    order_.clear(); coeff_.clear();
    NumMoments_ = order.size();
    if ( NumMoments_ != static_cast<uint>(coeff.size()) ) {
      coeff.resize(NumMoments_,one);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      order_.push_back((order[i] < two) ? two : order[i]);
      coeff_.push_back((coeff[i] < zero) ? one : coeff[i]);
    }
    // Build (approximate) positive function
    if ( list.get("Deviation Type","Upper") == "Upper" ) {
      positiveFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    }
    else {
      positiveFunction_ = Teuchos::rcp(new AbsoluteValue<Real>(list));
    }
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    if ( firstReset_ ) {
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      dualVector3_ = (x0->dual()).clone();
      dualVector4_ = (x0->dual()).clone();
      firstReset_  = false;
    }
    dualVector1_->zero(); dualVector2_->zero();
    dualVector3_->zero(); dualVector4_->zero();
    value_storage_.clear();
    gradient_storage_.clear();
    gradvec_storage_.clear();
    hessvec_storage_.clear();
    weights_.clear();
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
  } 
  
  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    value_storage_.push_back(val);
    weights_.push_back(weight);    
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_, ev(0), zero(0);
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    val = zero;
    Real diff(0), pf0(0), var(0);
    for ( uint i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        val += coeff_[p] * std::pow(pf0,order_[p]) * weights_[i];
      }
    }
    sampler.sumAll(&val,&var,1);
    // Return mean plus deviation
    return ev + var;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    RiskMeasure<Real>::g_->axpy(weight,g);
    value_storage_.push_back(val);
    gradient_storage_.push_back(g.clone());
    typename std::vector<Teuchos::RCP<Vector<Real> > >::iterator it = gradient_storage_.end();
    it--;
    (*it)->set(g);
    weights_.push_back(weight);    
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_, ev(0), zero(0), one(1);
    sampler.sumAll(&val,&ev,1);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector3_);
    // Compute deviation
    Real diff(0), pf0(0), pf1(0), c(0), ec(0), ecs(0);
    for ( uint i = 0; i < weights_.size(); i++ ) {
      c    = zero;
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        c += coeff_[p]*order_[p]*std::pow(pf0,order_[p]-one)*pf1;
      }
      ec += weights_[i]*c;
      dualVector1_->axpy(weights_[i]*c,*(gradient_storage_[i]));
    }
    sampler.sumAll(&ec,&ecs,1);
    dualVector3_->scale(one-ecs);
    sampler.sumAll(*dualVector1_,*dualVector2_);
    dualVector3_->plus(*dualVector2_);
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector3_);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    RiskMeasure<Real>::gv_  += weight * gv;
    RiskMeasure<Real>::g_->axpy(weight,g);
    RiskMeasure<Real>::hv_->axpy(weight,hv);
    value_storage_.push_back(val);
    gradient_storage_.push_back(g.clone());
    typename std::vector<Teuchos::RCP<Vector<Real> > >::iterator it = gradient_storage_.end();
    it--;
    (*it)->set(g);
    gradvec_storage_.push_back(gv);
    hessvec_storage_.push_back(hv.clone());
    it = hessvec_storage_.end();
    it--;
    (*it)->set(hv);
    weights_.push_back(weight);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    hv.zero();
    // Compute expected value
    std::vector<Real> myval(2), val(2);
    myval[0] = RiskMeasure<Real>::val_;
    myval[1] = RiskMeasure<Real>::gv_;
    sampler.sumAll(&myval[0],&val[0],2);
    Real ev = myval[0], egv = myval[1];
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector3_);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector4_);
    // Compute deviation
    Real diff(0), pf0(0), pf1(0), pf2(0), zero(0), one(1), two(2);
    Real cg(0), ecg(0), ecgs(0), ch(0), ech(0), echs(0);
    for ( uint i = 0; i < weights_.size(); i++ ) {
      cg   = zero;
      ch   = zero;
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        cg += coeff_[p]*order_[p]*(gradvec_storage_[i]-egv)*
                ((order_[p]-one)*std::pow(pf0,order_[p]-two)*pf1*pf1+
                std::pow(pf0,order_[p]-one)*pf2);
        ch += coeff_[p]*order_[p]*std::pow(pf0,order_[p]-one)*pf1;
      }
      ecg += weights_[i]*cg;
      ech += weights_[i]*ch;
      dualVector1_->axpy(weights_[i]*cg,*(gradient_storage_[i]));
      dualVector1_->axpy(weights_[i]*ch,*(hessvec_storage_[i]));
    }
    sampler.sumAll(&ech,&echs,1);
    dualVector4_->scale(one-echs);
    sampler.sumAll(&ecg,&ecgs,1);
    dualVector4_->axpy(-ecgs,*dualVector3_);
    sampler.sumAll(*dualVector1_,*dualVector2_);
    dualVector4_->plus(*dualVector2_);
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector4_);
  }
};

}

#endif
