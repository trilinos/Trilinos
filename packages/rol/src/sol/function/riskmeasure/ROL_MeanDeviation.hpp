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
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace ROL {

template<class Real>
class MeanDeviation : public RiskMeasure<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;

  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

  std::vector<Real> weights_;
  std::vector<Real> value_storage_;
  std::vector<Real> gradvec_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > gradient_storage_;
  std::vector<Teuchos::RCP<Vector<Real> > > hessvec_storage_;

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

  bool firstReset_;

  void initialize(void) {
    dev0_.clear(); dev1_.clear(); dev2_.clear(); dev3_.clear();
    des0_.clear(); des1_.clear(); des2_.clear(); des3_.clear();
    devp_.clear();
    gvp1_.clear(); gvp2_.clear(); gvp3_.clear();
    gvs1_.clear(); gvs2_.clear(); gvs3_.clear();

    dev0_.resize(NumMoments_,0.0); dev1_.resize(NumMoments_,0.0);
    dev2_.resize(NumMoments_,0.0); dev3_.resize(NumMoments_,0.0);
    des0_.resize(NumMoments_,0.0); des1_.resize(NumMoments_,0.0);
    des2_.resize(NumMoments_,0.0); des3_.resize(NumMoments_,0.0);
    devp_.resize(NumMoments_,0.0);
    gvp1_.resize(NumMoments_,0.0); gvp2_.resize(NumMoments_,0.0);
    gvp3_.resize(NumMoments_,0.0);
    gvs1_.resize(NumMoments_,0.0); gvs2_.resize(NumMoments_,0.0);
    gvs3_.resize(NumMoments_,0.0);
  }

  void clear(void) {
    dev0_.assign(NumMoments_,0.0); dev1_.assign(NumMoments_,0.0);
    dev2_.assign(NumMoments_,0.0); dev3_.assign(NumMoments_,0.0);
    des0_.assign(NumMoments_,0.0); des1_.assign(NumMoments_,0.0);
    des2_.assign(NumMoments_,0.0); des3_.assign(NumMoments_,0.0);
    devp_.assign(NumMoments_,0.0);
    gvp1_.assign(NumMoments_,0.0); gvp2_.assign(NumMoments_,0.0);
    gvp3_.assign(NumMoments_,0.0);
    gvs1_.assign(NumMoments_,0.0); gvs2_.assign(NumMoments_,0.0);
    gvs3_.assign(NumMoments_,0.0);

    value_storage_.clear();
    gradient_storage_.clear();
    gradvec_storage_.clear();
    hessvec_storage_.clear();
    weights_.clear();
  }

public:

  MeanDeviation( Real order, Real coeff,
                 Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    order_.clear(); coeff_.clear();
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
    NumMoments_ = order_.size();
    initialize();
  }

  MeanDeviation( std::vector<Real> &order, std::vector<Real> &coeff, 
                 Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    order_.clear(); coeff_.clear();
    NumMoments_ = order.size();
    if ( NumMoments_ != coeff.size() ) {
      coeff.resize(NumMoments_,1.0);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      order_.push_back((order[i] < 2.0) ? 2.0 : order[i]);
      coeff_.push_back((coeff[i] < 0.0) ? 1.0 : coeff[i]);
    }
    initialize();
  }

  MeanDeviation( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Deviation");
    // Get data from parameter list
    Teuchos::Array<Real> order
      = Teuchos::getArrayFromStringParameter<double>(list,"Orders");
    Teuchos::Array<Real> coeff
      = Teuchos::getArrayFromStringParameter<double>(list,"Coefficients");
    // Check inputs
    NumMoments_ = order.size();
    order_.clear(); coeff_.clear();
    if ( NumMoments_ != static_cast<uint>(coeff.size()) ) {
      coeff.resize(NumMoments_,1.0);
    }
    for ( uint i = 0; i < NumMoments_; i++ ) {
      order_.push_back((order[i] < 2.0) ? 2.0 : order[i]);
      coeff_.push_back((coeff[i] < 0.0) ? 1.0 : coeff[i]);
    }
    // Build (approximate) positive function
    if ( list.get("Deviation Type","Upper") == "Upper" ) {
      positiveFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    }
    else {
      positiveFunction_ = Teuchos::rcp(new AbsoluteValue<Real>(list));
    }
    initialize();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    if ( firstReset_ ) {
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      firstReset_  = false;
    }
    dualVector1_->zero(); dualVector2_->zero();
    clear();
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
    Real val = RiskMeasure<Real>::val_, ev = 0.0;
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, dev = 0.0;
    for ( uint i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += std::pow(pf0,order_[p]) * weights_[i];
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev += coeff_[p]*std::pow(des0_[p],1.0/order_[p]);
    }
    // Return mean plus deviation
    return ev + dev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_, ev = 0.0;
    sampler.sumAll(&val,&ev,1);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 0.0;
    for ( uint i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += weights_[i] * std::pow(pf0,order_[p]);
        dev1_[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf1;
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    sampler.sumAll(&dev1_[0],&des1_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      dev0_[p] = std::pow(des0_[p],1.0-1.0/order_[p]);
    }
    // Compute derivative
    for ( uint i = 0; i < weights_.size(); i++ ) {
      c    = 0.0;
      diff = value_storage_[i]-ev;
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        if ( dev0_[p] > 0.0 ) {
          c += coeff_[p]/dev0_[p] * (std::pow(pf0,order_[p]-1.0)*pf1 - des1_[p]);
        }
      }
      dualVector1_->axpy(weights_[i]*c,*(gradient_storage_[i]));
    }
    dualVector1_->axpy(1.0,*(RiskMeasure<Real>::g_));
    sampler.sumAll(*dualVector1_,*dualVector2_);
    // Set RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector2_);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    // Compute expected value
    Real val = RiskMeasure<Real>::val_, ev = 0.0;
    sampler.sumAll(&val,&ev,1);
    Real gv  = RiskMeasure<Real>::gv_, egv = 0.0;
    sampler.sumAll(&gv,&egv,1);
    // Compute deviation
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0;
    Real cg = 0.0, ch = 0.0, diff1 = 0.0, diff2 = 0.0, diff3 = 0.0;
    for ( uint i = 0; i < weights_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        dev0_[p] += weights_[i] * std::pow(pf0,order_[p]);
        dev1_[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf1;
        dev2_[p] += weights_[i] * std::pow(pf0,order_[p]-2.0) * pf1 * pf1;
        dev3_[p] += weights_[i] * std::pow(pf0,order_[p]-1.0) * pf2;
      }
    }
    sampler.sumAll(&dev0_[0],&des0_[0],NumMoments_);
    sampler.sumAll(&dev1_[0],&des1_[0],NumMoments_);
    sampler.sumAll(&dev2_[0],&des2_[0],NumMoments_);
    sampler.sumAll(&dev3_[0],&des3_[0],NumMoments_);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      devp_[p] = std::pow(des0_[p],2.0-1.0/order_[p]);
      dev0_[p] = std::pow(des0_[p],1.0-1.0/order_[p]);
    }
    for ( uint i = 0; i < value_storage_.size(); i++ ) {
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        gvp1_[p] += weights_[i] * (std::pow(pf0,order_[p]-1.0)*pf1-des1_[p]) *
                     (gradvec_storage_[i] - egv);
        gvp2_[p] += weights_[i] * (std::pow(pf0,order_[p]-2.0)*pf1*pf1-des2_[p]) *
                     (gradvec_storage_[i] - egv);
        gvp3_[p] += weights_[i] * (std::pow(pf0,order_[p]-1.0)*pf2-des3_[p]) *
                     (gradvec_storage_[i] - egv);
      }
    }
    sampler.sumAll(&gvp1_[0],&gvs1_[0],NumMoments_);
    sampler.sumAll(&gvp2_[0],&gvs2_[0],NumMoments_);
    sampler.sumAll(&gvp3_[0],&gvs3_[0],NumMoments_);
    // Compute derivative
    for ( uint i = 0; i < weights_.size(); i++ ) {
      cg   = 1.0;
      ch   = 0.0;
      diff = value_storage_[i]-ev;
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      for ( uint p = 0; p < NumMoments_; p++ ) {
        if ( dev0_[p] > 0.0 ) {
          diff1 = std::pow(pf0,order_[p]-1.0)*pf1-des1_[p];
          diff2 = std::pow(pf0,order_[p]-2.0)*pf1*pf1*(gradvec_storage_[i]-egv)-gvs2_[p];
          diff3 = std::pow(pf0,order_[p]-1.0)*pf2*(gradvec_storage_[i]-egv)-gvs3_[p];
          cg   += coeff_[p]*diff1/dev0_[p];
          ch   += coeff_[p]*(((order_[p]-1.0)*diff2+diff3)/dev0_[p] -
                    (order_[p]-1.0)*gvs1_[p]*diff1/devp_[p]);
        }
      }
      dualVector1_->axpy(weights_[i]*ch,*(gradient_storage_[i]));
      dualVector1_->axpy(weights_[i]*cg,*(hessvec_storage_[i]));
    }
    sampler.sumAll(*dualVector1_,*dualVector2_);
    // Fill RiskVector
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector2_);
  }
};

}

#endif
