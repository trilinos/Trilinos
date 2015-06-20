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

namespace ROL {

template<class Real>
class MeanVariance : public RiskMeasure<Real> {
private:
  std::vector<Real> order_;
  std::vector<Real> coeff_;

  std::vector<Real> weights_;
  std::vector<Real> value_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > gradient_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > hessvec_storage_;
  std::vector<Real> gradvec_storage_;

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

public:
  MeanVariance( Real order, Real coeff,
                Teuchos::RCP<PositiveFunction<Real> > &pf ) : positiveFunction_(pf) {
    order_.clear();
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.clear();
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
  }
  MeanVariance( std::vector<Real> &order, std::vector<Real> &coeff, 
                Teuchos::RCP<PositiveFunction<Real> > &pf ) : positiveFunction_(pf) {
    order_.clear();
    coeff_.clear();
    if ( order.size() != coeff.size() ) {
      coeff.resize(order.size(),1.0);
    }
    for ( unsigned i = 0; i < order.size(); i++ ) {
      order_.push_back((order[i] < 2.0) ? 2.0 : order[i]);
      coeff_.push_back((coeff[i] < 0.0) ? 1.0 : coeff[i]);
    }
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::gv_  = 0.0;
    RiskMeasure<Real>::g_   = x.clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_  = x.clone(); RiskMeasure<Real>::hv_->zero();
    (this->value_storage_).clear();
    (this->gradient_storage_).clear();
    (this->gradvec_storage_).clear();
    (this->hessvec_storage_).clear();
    (this->weights_).clear();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false); 
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::gv_  = 0.0;
    RiskMeasure<Real>::g_   = x.clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_  = x.clone(); RiskMeasure<Real>::hv_->zero();
    (this->value_storage_).clear();
    (this->gradient_storage_).clear();
    (this->gradvec_storage_).clear();
    (this->hessvec_storage_).clear();
    (this->weights_).clear();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false);
    v0 = Teuchos::rcp(&const_cast<Vector<Real> &>(v),false);
  } 
  
  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    this->value_storage_.push_back(val);
    this->weights_.push_back(weight);    
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    RiskMeasure<Real>::g_->axpy(weight,g);
    this->value_storage_.push_back(val);
    this->gradient_storage_.push_back(g.clone());
    typename std::vector<Teuchos::RCP<Vector<Real> > >::iterator it = (this->gradient_storage_).end();
    it--;
    (*it)->set(g);
    this->weights_.push_back(weight);    
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    RiskMeasure<Real>::gv_  += weight * gv;
    RiskMeasure<Real>::g_->axpy(weight,g);
    RiskMeasure<Real>::hv_->axpy(weight,hv);
    this->value_storage_.push_back(val);
    this->gradient_storage_.push_back(g.clone());
    typename std::vector<Teuchos::RCP<Vector<Real> > >::iterator it = (this->gradient_storage_).end();
    it--;
    (*it)->set(g);
    this->gradvec_storage_.push_back(gv);
    this->hessvec_storage_.push_back(hv.clone());
    it = (this->hessvec_storage_).end();
    it--;
    (*it)->set(hv);
    this->weights_.push_back(weight);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    val = 0.0;
    Real diff = 0.0, pf0 = 0.0, var = 0.0;
    for ( unsigned i = 0; i < this->weights_.size(); i++ ) {
      diff = this->value_storage_[i]-ev;
      pf0  = this->positiveFunction_->evaluate(diff,0);
      for ( unsigned p = 0; p < this->order_.size(); p++ ) {
        val += this->coeff_[p] * std::pow(pf0,this->order_[p]) * this->weights_[i];
      }
    }
    sampler.sumAll(&val,&var,1);
    // Return mean plus deviation
    return ev + var;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    g.zero();
    // Compute expected value
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    sampler.sumAll(*(RiskMeasure<Real>::g_),g);
    // Compute deviation
    Teuchos::RCP<Vector<Real> > gs   = g.clone(); gs->zero();
    Teuchos::RCP<Vector<Real> > gtmp = g.clone(); gtmp->zero();
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 0.0, ec = 0.0, ecs = 0.0;
    for ( unsigned i = 0; i < this->weights_.size(); i++ ) {
      c    = 0.0;
      diff = this->value_storage_[i]-ev;
      pf0  = this->positiveFunction_->evaluate(diff,0);
      pf1  = this->positiveFunction_->evaluate(diff,1);
      for ( unsigned p = 0; p < this->order_.size(); p++ ) {
        c += this->coeff_[p]*this->order_[p]*std::pow(pf0,this->order_[p]-1.0)*pf1;
      }
      ec += this->weights_[i]*c;
      gtmp->axpy(this->weights_[i]*c,*(this->gradient_storage_[i]));
    }
    sampler.sumAll(&ec,&ecs,1);
    g.scale(1.0-ecs);
    sampler.sumAll(*gtmp,*gs);
    g.plus(*gs);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    hv.zero();
    // Compute expected value
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    Real gv  = RiskMeasure<Real>::gv_;
    Real egv = 0.0;
    sampler.sumAll(&gv,&egv,1);
    Teuchos::RCP<Vector<Real> > g = hv.clone();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*g);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),hv);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0;  
    Real cg = 0.0, ecg = 0.0, ecgs = 0.0, ch = 0.0, ech = 0.0, echs = 0.0;
    Teuchos::RCP<Vector<Real> > htmp = hv.clone(); htmp->zero();
    Teuchos::RCP<Vector<Real> > hs   = hv.clone(); hs->zero();
    for ( unsigned i = 0; i < this->weights_.size(); i++ ) {
      cg   = 0.0;
      ch   = 0.0;
      diff = this->value_storage_[i]-ev;
      pf0  = this->positiveFunction_->evaluate(diff,0);
      pf1  = this->positiveFunction_->evaluate(diff,1);
      pf2  = this->positiveFunction_->evaluate(diff,2);
      for ( unsigned p = 0; p < this->order_.size(); p++ ) {
        cg += this->coeff_[p]*this->order_[p]*(this->gradvec_storage_[i]-egv)*
                ((this->order_[p]-1.0)*std::pow(pf0,this->order_[p]-2.0)*pf1*pf1+
                std::pow(pf0,this->order_[p]-1.0)*pf2);
        ch += this->coeff_[p]*this->order_[p]*std::pow(pf0,this->order_[p]-1.0)*pf1;
      }
      ecg += this->weights_[i]*cg;
      ech += this->weights_[i]*ch;
      htmp->axpy(this->weights_[i]*cg,*(this->gradient_storage_[i]));
      htmp->axpy(this->weights_[i]*ch,*(this->hessvec_storage_[i]));
    }
    sampler.sumAll(&ech,&echs,1);
    hv.scale(1.0-echs);
    sampler.sumAll(&ecg,&ecgs,1);
    hv.axpy(-ecgs,*g);
    sampler.sumAll(*htmp,*hs);
    hv.plus(*hs);
  }
};

}

#endif
