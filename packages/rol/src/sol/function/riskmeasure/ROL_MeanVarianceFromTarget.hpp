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
#include "ROL_PlusFunction.hpp"
#include "ROL_AbsoluteValue.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

namespace ROL {

template<class Real>
class MeanVarianceFromTarget : public RiskMeasure<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:

  Teuchos::RCP<PositiveFunction<Real> > positiveFunction_;

  std::vector<Real> target_;
  std::vector<Real> order_;
  std::vector<Real> coeff_;
  uint NumMoments_;

public:

  MeanVarianceFromTarget( Real target, Real order, Real coeff,
                          Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
    Real zero(0), one(1), two(2);
    target_.clear(); order_.clear(); coeff_.clear();
    target_.push_back(target);
    order_.push_back((order < two) ? two : order);
    coeff_.push_back((coeff < zero) ? one : coeff);
    NumMoments_ = order_.size();
  }

  MeanVarianceFromTarget( std::vector<Real> &target, std::vector<Real> &order, std::vector<Real> &coeff, 
                          Teuchos::RCP<PositiveFunction<Real> > &pf )
    : RiskMeasure<Real>(), positiveFunction_(pf) {
    Real zero(0), one(1), two(2);
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
      coeff_.push_back((coeff[i] < zero) ? one : coeff[i]);
    }
  }
  
  MeanVarianceFromTarget( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>() {
    Real zero(0), one(1), two(2);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean Plus Variance From Target");
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
    // Build (approximate) positive function
    if ( list.get("Deviation Type","Upper") == "Upper" ) {
      positiveFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    }
    else {
      positiveFunction_ = Teuchos::rcp(new AbsoluteValue<Real>(list));
    }
  }
  
  void update(const Real val, const Real weight) {
    Real diff(0), pf0(0);
    RiskMeasure<Real>::val_ += weight * val;
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      RiskMeasure<Real>::val_ += weight * coeff_[p] * std::pow(pf0,order_[p]);
    }
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real diff(0), pf0(0), pf1(0), c(1), one(1);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val-target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      c   += order_[p]*coeff_[p]*std::pow(pf0,order_[p]-one)*pf1;
    }
    (RiskMeasure<Real>::g_)->axpy(weight * c,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real diff(0), pf0(0), pf1(0), pf2(0), p1(0), p2(0), ch(1), cg(0), one(1), two(2);
    for ( uint p = 0; p < NumMoments_; p++ ) {
      diff = val - target_[p];
      pf0  = positiveFunction_->evaluate(diff,0);
      pf1  = positiveFunction_->evaluate(diff,1);
      pf2  = positiveFunction_->evaluate(diff,2);
      //p0   = std::pow(pf0,order_[p]);
      p1   = std::pow(pf0,order_[p]-one);
      p2   = std::pow(pf0,order_[p]-two);
      cg  += order_[p]*coeff_[p]*gv*( (order_[p]-one)*p2*pf1*pf1 + p1*pf2 );
      ch  += order_[p]*coeff_[p]*p1*pf1;
    }
    RiskMeasure<Real>::hv_->axpy(weight*cg,g);
    RiskMeasure<Real>::hv_->axpy(weight*ch,hv);
  }
};

}

#endif
