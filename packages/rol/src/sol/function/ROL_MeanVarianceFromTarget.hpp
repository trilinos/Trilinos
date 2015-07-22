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

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;

public:

  MeanVarianceFromTarget( Real target, Real order, Real coeff,
                          Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
    target_.clear(); order_.clear(); coeff_.clear();
    target_.push_back(target);
    order_.push_back((order < 2.0) ? 2.0 : order);
    coeff_.push_back((coeff < 0.0) ? 1.0 : coeff);
  }

  MeanVarianceFromTarget( std::vector<Real> &target, std::vector<Real> &order, std::vector<Real> &coeff, 
                          Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
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
  }
  
  void update(const Real val, const Real weight) {
    Real diff = 0.0, pf0 = 0.0;
    RiskMeasure<Real>::val_ += weight * val;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      RiskMeasure<Real>::val_ += weight * coeff_[p] * std::pow(pf0,order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, c = 1.0;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      c   += order_[p]*coeff_[p]*std::pow(pf0,order_[p]-1.0)*pf1;
    }
    (RiskMeasure<Real>::g_)->axpy(weight * c,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff = 0.0, pf0 = 0.0, pf1 = 0.0, pf2 = 0.0, p1 = 0.0, p2 = 0.0, ch = 1.0, cg = 0.0;
    for ( unsigned p = 0; p < order_.size(); p++ ) {
      diff = val - target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      //p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-1.0);
      p2   = std::pow(pf0,order_[p]-2.0);
      cg  += order_[p]*coeff_[p]*gv*( (order_[p]-1.0)*p2*pf1*pf1 + p1*pf2 );
      ch  += order_[p]*coeff_[p]*p1*pf1;
    }
    RiskMeasure<Real>::hv_->axpy(weight*cg,g);
    RiskMeasure<Real>::hv_->axpy(weight*ch,hv);
  }
};

}

#endif
