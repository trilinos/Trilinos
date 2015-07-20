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

#ifndef ROL_MEANDEVIATION_HPP
#define ROL_MEANDEVIATION_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PositiveFunction.hpp"

namespace ROL {

template<class Real>
class MeanDeviation : public RiskMeasure<Real> {
private:

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  std::vector<Real> order_;
  std::vector<Real> coeff_;

  std::vector<Real> weights_;
  std::vector<Real> value_storage_;
  std::vector<Real> gradvec_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > gradient_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > hessvec_storage_;

public:

  MeanDeviation( Real order, Real coeff,
                 Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
    order_.clear(); coeff_.clear();
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
  }

  MeanDeviation( std::vector<Real> &order, std::vector<Real> &coeff, 
                 Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
    order_.clear(); coeff_.clear();
    if ( order.size() != coeff.size() ) {
      coeff.resize(order.size(),1.0);
    }
    for ( unsigned i = 0; i < order.size(); i++ ) {
      order_.push_back((order[i] < 2.0) ? 2.0 : order[i]);
      coeff_.push_back((coeff[i] < 0.0) ? 1.0 : coeff[i]);
    }
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    value_storage_.clear();
    gradient_storage_.clear();
    gradvec_storage_.clear();
    hessvec_storage_.clear();
    weights_.clear();
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    RiskMeasure<Real>::reset(x0,x,v0,v);
    value_storage_.clear();
    gradient_storage_.clear();
    gradvec_storage_.clear();
    hessvec_storage_.clear();
    weights_.clear();
  }
  
  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * val;
    value_storage_.push_back(val);
    weights_.push_back(weight);    
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

  Real getValue(SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, dev = 0.0;
    std::vector<Real> devp(order_.size(),0.0);
    std::vector<Real> devs(order_.size(),0.0);
    for ( unsigned i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        devp[p] += std::pow(pf0,order_[p]) * weights_[i];
      }
    }
    sampler.sumAll(&devp[0],&devs[0],devp.size());
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      dev += coeff_[p]*std::pow(devs[p],1.0/order_[p]);
    }
    // Return mean plus deviation
    return ev + dev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    g.zero();
    // Compute expected value
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 0.0;
    std::vector<Real> dev0(order_.size(),0.0);
    std::vector<Real> des0(order_.size(),0.0);
    std::vector<Real> dev1(order_.size(),0.0);
    std::vector<Real> des1(order_.size(),0.0);
    for ( unsigned i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        dev0[p] += weights_[i] * std::pow(pf0,order_[p]);
        dev1[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf1;
      }
    }
    sampler.sumAll(&dev0[0],&des0[0],dev0.size());
    sampler.sumAll(&dev1[0],&des1[0],dev1.size());
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      dev0[p] = std::pow(des0[p],1.0-1.0/order_[p]);
    }
    // Compute derivative
    Teuchos::RCP<Vector<Real> > gtmp = g.clone();
    gtmp->zero();
    for ( unsigned i = 0; i < weights_.size(); i++ ) {
      c    = 0.0;
      diff = value_storage_[i]-ev;
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        if ( dev0[p] > 0.0 ) {
          c += coeff_[p]/dev0[p] * (std::pow(pf0,order_[p]-1.0)*pf1 - des1[p]);
        }
      }
      gtmp->axpy(weights_[i]*c,*(gradient_storage_[i]));
    }
    gtmp->axpy(1.0,*(RiskMeasure<Real>::g_));
    sampler.sumAll(*gtmp,g);
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
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0;
    Real cg = 0.0, ch = 0.0, diff1 = 0.0, diff2 = 0.0, diff3 = 0.0;
    std::vector<Real> dev0(order_.size(),0.0);
    std::vector<Real> dev1(order_.size(),0.0);
    std::vector<Real> dev2(order_.size(),0.0);
    std::vector<Real> dev3(order_.size(),0.0);
    std::vector<Real> des0(order_.size(),0.0);
    std::vector<Real> des1(order_.size(),0.0);
    std::vector<Real> des2(order_.size(),0.0);
    std::vector<Real> des3(order_.size(),0.0);
    std::vector<Real> devp(order_.size(),0.0);
    std::vector<Real> gvp1(order_.size(),0.0);
    std::vector<Real> gvp2(order_.size(),0.0);
    std::vector<Real> gvp3(order_.size(),0.0);
    std::vector<Real> gvs1(order_.size(),0.0);
    std::vector<Real> gvs2(order_.size(),0.0);
    std::vector<Real> gvs3(order_.size(),0.0);
    for ( unsigned i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        dev0[p] += weights_[i] * std::pow(pf0,order_[p]);
        dev1[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf1;
        dev2[p] += weights_[i] * std::pow(pf0,order_[p]-2.0) * pf1 * pf1;
        dev3[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf2;
      }
    }
    sampler.sumAll(&dev0[0],&des0[0],dev0.size());
    sampler.sumAll(&dev1[0],&des1[0],dev1.size());
    sampler.sumAll(&dev2[0],&des2[0],dev2.size());
    sampler.sumAll(&dev3[0],&des3[0],dev3.size());
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      devp[p] = std::pow(des0[p],2.0-1.0/order_[p]);
      dev0[p] = std::pow(des0[p],1.0-1.0/order_[p]);
    }
    for ( unsigned i = 0; i < value_storage_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        gvp1[p] += weights_[i] * (std::pow(pf0,order_[p]-1.0)*pf1-des1[p]) *
                     (gradvec_storage_[i] - egv);
        gvp2[p] += weights_[i] * (std::pow(pf0,order_[p]-2.0)*pf1*pf1-des2[p]) *
                     (gradvec_storage_[i] - egv);
        gvp3[p] += weights_[i] * (std::pow(pf0,order_[p]-1.0)*pf2-des3[p]) *
                     (gradvec_storage_[i] - egv);
      }
    }
    sampler.sumAll(&gvp1[0],&gvs1[0],gvp1.size());
    sampler.sumAll(&gvp2[0],&gvs2[0],gvp2.size());
    sampler.sumAll(&gvp3[0],&gvs3[0],gvp3.size());
    // Compute derivative
    Teuchos::RCP<Vector<Real> > htmp = hv.clone();
    htmp->zero();
    for ( unsigned i = 0; i < weights_.size(); i++ ) {
      cg   = 1.0;
      ch   = 0.0;
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        if ( dev0[p] > 0.0 ) {
          diff1 = std::pow(pf0,order_[p]-1.0)*pf1-des1[p];
          diff2 = std::pow(pf0,order_[p]-2.0)*pf1*pf1*(gradvec_storage_[i]-egv)-gvs2[p];
          diff3 = std::pow(pf0,order_[p]-1.0)*pf2*(gradvec_storage_[i]-egv)-gvs3[p];
          cg   += coeff_[p]*diff1/dev0[p];
          ch   += coeff_[p]*(((order_[p]-1.0)*diff2+diff3)/dev0[p] -
                    (order_[p]-1.0)*gvs1[p]*diff1/devp[p]);
        }
      }
      htmp->axpy(weights_[i]*ch,*(gradient_storage_[i]));
      htmp->axpy(weights_[i]*cg,*(hessvec_storage_[i]));
    }
    sampler.sumAll(*htmp,hv);
  }
};

}

#endif
