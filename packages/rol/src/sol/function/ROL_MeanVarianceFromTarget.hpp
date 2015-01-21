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

#ifndef ROL_MEANVARIANCEFROMTARGET_HPP
#define ROL_MEANVARIANCEFROMTARGET_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PositiveFunction.hpp"

namespace ROL {

template<class Real>
class MeanVarianceFromTarget : public RiskMeasure<Real> {
private:
  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

public:
  MeanVarianceFromTarget( Real target, Real order, Real coeff,
                          Teuchos::RCP<PositiveFunction<Real> > &pf ) : positiveFunction_(pf) {
    target_.clear();
    target_.push_back(target);
    order_.clear();
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.clear();
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
  }
  MeanVarianceFromTarget( std::vector<Real> &target, std::vector<Real> &order, std::vector<Real> &coeff, 
                          Teuchos::RCP<PositiveFunction<Real> > &pf ) : positiveFunction_(pf) {
    target_.clear();
    order_.clear();
    coeff_.clear();
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
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::gv_  = 0.0;
    RiskMeasure<Real>::g_  = x.clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x.clone(); RiskMeasure<Real>::hv_->zero();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false); 
  }
    
  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    RiskMeasure<Real>::val_ = 0.0;
    RiskMeasure<Real>::gv_  = 0.0;
    RiskMeasure<Real>::g_  = x.clone(); RiskMeasure<Real>::g_->zero();
    RiskMeasure<Real>::hv_ = x.clone(); RiskMeasure<Real>::hv_->zero();
    x0 = Teuchos::rcp(&const_cast<Vector<Real> &>(x),false);
    v0 = Teuchos::rcp(&const_cast<Vector<Real> &>(v),false);
  } 
  
  void update(const Real val, const Real weight) {
    Real diff = 0.0, pf0 = 0.0;
    RiskMeasure<Real>::val_ += weight * val;
    for ( unsigned p = 0; p < this->order_.size(); p++ ) {
      diff = val-this->target_[p];
      pf0  = this->positiveFunction_->evaluate(diff,0);
      RiskMeasure<Real>::val_ += weight * this->coeff_[p] * std::pow(pf0,this->order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 1.0;
    for ( unsigned p = 0; p < this->order_.size(); p++ ) {
      diff = val-this->target_[p];
      pf0  = this->positiveFunction_->evaluate(diff,0);
      pf1  = this->positiveFunction_->evaluate(diff,1);
      c   += this->order_[p]*this->coeff_[p]*std::pow(pf0,this->order_[p]-1.0)*pf1;
    }
    (RiskMeasure<Real>::g_)->axpy(weight * c,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0, p0 = 0.0, p1 = 0.0, p2 = 0.0, ch = 1.0, cg = 0.0;
    for ( unsigned p = 0; p < this->order_.size(); p++ ) {
      diff = val - this->target_[p];
      pf0  = this->positiveFunction_->evaluate(diff,0);
      pf1  = this->positiveFunction_->evaluate(diff,1);
      pf2  = this->positiveFunction_->evaluate(diff,2);
      p0   = std::pow(pf0,this->order_[p]);
      p1   = std::pow(pf0,this->order_[p]-1.0);
      p2   = std::pow(pf0,this->order_[p]-2.0);
      cg  += this->order_[p]*this->coeff_[p]*gv*( (this->order_[p]-1.0)*p2*pf1*pf1 + p1*pf2 );
      ch  += this->order_[p]*this->coeff_[p]*p1*pf1;
    }
    RiskMeasure<Real>::hv_->axpy(weight*cg,g);
    RiskMeasure<Real>::hv_->axpy(weight*ch,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_;
    Real ev  = 0.0;
    sampler.sumAll(&val,&ev,1);
    return ev;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(RiskMeasure<Real>::g_),g);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    sampler.sumAll(*(RiskMeasure<Real>::hv_),hv);
  }
};

}

#endif
