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

namespace ROL {

template<class Real>
class MeanDeviationFromTarget : public RiskMeasure<Real> {
private:
  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  std::vector<Real> pval_; 
  std::vector<Real> pgv_; 

  std::vector<Teuchos::RCP<Vector<Real> > > pg0_;
  std::vector<Teuchos::RCP<Vector<Real> > > pg_;
  std::vector<Teuchos::RCP<Vector<Real> > > phv_;

  bool firstReset_;

public:

  MeanDeviationFromTarget( Real target, Real order, Real coeff,
                           Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    // Initialize storage for problem data
    target_.clear(); order_.clear(); coeff_.clear();
    target_.push_back(target);
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
    // Initialize additional storage
    pg_.clear(); pg0_.clear(); phv_.clear(); pval_.clear(); pgv_.clear();
    pg_.resize(order_.size());
    pg0_.resize(order_.size());
    phv_.resize(order_.size());
    pval_.resize(order_.size(),0.0);
    pgv_.resize(order_.size(),0.0);
  }

  MeanDeviationFromTarget( std::vector<Real> &target, std::vector<Real> &order, std::vector<Real> &coeff, 
                           Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf), firstReset_(true) {
    // Initialize storage for problem data
    target_.clear(); order_.clear(); coeff_.clear();
    if ( order.size() != target.size() ) {
      target.resize(order.size(),0.0);
    }
    if ( order.size() != coeff.size() ) {
      coeff.resize(order.size(),1.0);
    }
    for ( unsigned i = 0; i < order.size(); i++ ) {
      target_.push_back(target[i]);
      order_.push_back((order[i] < 2.0) ? 2.0 : order[i]);
      coeff_.push_back((coeff[i] < 0.0) ? 1.0 : coeff[i]);
    }
    // Initialize additional storage
    pg_.clear(); pg0_.clear(); phv_.clear(); pval_.clear(); pgv_.clear();
    pg_.resize(order_.size());
    pg0_.resize(order_.size());
    phv_.resize(order_.size());
    pval_.resize(order_.size(),0.0);
    pgv_.resize(order_.size(),0.0);
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    if (firstReset_) {
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        pg0_[p] = (x.dual()).clone();
        pg_[p]  = (x.dual()).clone();
        phv_[p] = (x.dual()).clone();
      }
      firstReset_ = false;
    }
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      pg0_[p]->zero(); pg_[p]->zero(); phv_[p]->zero();
      pval_[p] = 0.0; pgv_[p] = 0.0;
    }
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    RiskMeasure<Real>::reset(x0,x,v0,v);
    if (firstReset_) {
      for ( unsigned p = 0; p < order_.size(); p++ ) {
        pg0_[p] = (x.dual()).clone();
        pg_[p]  = (x.dual()).clone();
        phv_[p] = (x.dual()).clone();
      }
      firstReset_ = false;
    }
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      pg0_[p]->zero(); pg_[p]->zero(); phv_[p]->zero();
      pval_[p] = 0.0; pgv_[p] = 0.0;
    }
  }
  
  void update(const Real val, const Real weight) {
    Real diff = 0.0, pf0 = 0.0;
    RiskMeasure<Real>::val_ += weight * val;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 0.0;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val-target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      c    = std::pow(pf0,order_[p]-1.0) * pf1;
      (pg_[p])->axpy(weight * c,g);
      pval_[p] += weight * std::pow(pf0,order_[p]);
    }
    RiskMeasure<Real>::g_->axpy(weight,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0, p0 = 0.0, p1 = 0.0, p2 = 0.0, c = 0.0;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val - target_[p];
      pf0 = positiveFunction_->evaluate(diff,0);
      pf1 = positiveFunction_->evaluate(diff,1);
      pf2 = positiveFunction_->evaluate(diff,2);
      p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-1.0);
      p2   = std::pow(pf0,order_[p]-2.0);
      c    = -(order_[p]-1.0)*p1*pf1;
      pg0_[p]->axpy(weight*c,g);
      c    = gv*((order_[p]-1.0)*p2*pf1*pf1 + p1*pf2);
      pg_[p]->axpy(weight*c,g);
      c    = p1*pf1;
      phv_[p]->axpy(weight*c,hv);
      pval_[p] += weight*p0;
      pgv_[p]  += weight*p1*pf1*gv;
    }
    RiskMeasure<Real>::hv_->axpy(weight,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_;
    Real dev = 0.0;
    sampler.sumAll(&val,&dev,1);
    std::vector<Real> pval_sum(pval_.size());
    sampler.sumAll(&(pval_)[0],&pval_sum[0],pval_.size());
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      dev += coeff_[p] * std::pow(pval_sum[p],1.0/order_[p]);
    }
    return dev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(RiskMeasure<Real>::g_),g);
    std::vector<Real> pval_sum(pval_.size());
    sampler.sumAll(&(pval_)[0],&pval_sum[0],pval_.size());
    Teuchos::RCP<Vector<Real> > pg;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      if ( pval_sum[p] > 0.0 ) {
        pg = (pg_[p])->clone();
        sampler.sumAll(*(pg_[p]),*pg);
        g.axpy(coeff_[p]/std::pow(pval_sum[p],1.0-1.0/order_[p]),*pg);
      }
    }
  }
  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(RiskMeasure<Real>::hv_),hv);
    std::vector<Real> pval_sum(pval_.size());
    sampler.sumAll(&(pval_)[0],&pval_sum[0],pval_.size());
    std::vector<Real> pgv_sum(pgv_.size());
    sampler.sumAll(&(pgv_)[0],&pgv_sum[0],pgv_.size());
    Real c = 0.0;
    Teuchos::RCP<Vector<Real> > pg, pg0, phv;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      if ( pval_sum[p] > 0.0 ) {
        pg  = (pg_[p])->clone();
        sampler.sumAll(*(pg_[p]),*pg);
        pg0 = (pg0_[p])->clone();
        sampler.sumAll(*(pg0_[p]),*pg0);
        phv = (phv_[p])->clone();
        sampler.sumAll(*(phv_[p]),*phv);
        c = coeff_[p]*(pgv_sum[p]/std::pow(pval_sum[p],2.0-1.0/order_[p]));
        hv.axpy(c,*pg0);
        c = coeff_[p]/std::pow(pval_sum[p],1.0-1.0/order_[p]);
        hv.axpy(c,*pg);
        hv.axpy(c,*phv);
      }
    }
  }
};

}

#endif
