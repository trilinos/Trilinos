// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUNCATEDEXPONENTIAL_HPP
#define ROL_TRUNCATEDEXPONENTIAL_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class TruncatedExponential : public Distribution<Real> {
private:
  Real a_;
  Real b_;
  Real scale_;
  Real expa_;
  Real expb_;
  Real diff_;
  Real coeff_;

  size_t compute_coeff(const size_t m, const size_t k) const {
    if ( k == m || m == 0 || m == 1 ) {
      return 1;
    }
    size_t val = 1;
    for (size_t i = k; i < m; i++) {
      val *= (i+1);
    }
    return val;
  }

public: 
  TruncatedExponential(const Real a = 0., const Real b = 1., const Real scale = 1.)
    : a_(std::min(a,b)), b_(std::max(a,b)), scale_((scale>0.) ? scale : 1.) {
    expa_  = std::exp(-scale_*a_);
    expb_  = std::exp(-scale_*b_);
    diff_  = expa_ - expb_;
    coeff_ = scale_/diff_;
  }

  TruncatedExponential(ROL::ParameterList &parlist) {
    ROL::ParameterList TElist
      = parlist.sublist("SOL").sublist("Distribution").sublist("Truncated Exponential");
    a_ = TElist.get("Lower Bound",0.);
    b_ = TElist.get("Upper Bound",1.);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);
    scale_ = TElist.get("Scale",1.);
    scale_ = (scale_ > 0.) ? scale_ : 1.;
    expa_  = std::exp(-scale_*a_);
    expb_  = std::exp(-scale_*b_);
    diff_  = expa_ - expb_;
    coeff_ = scale_/diff_;
  }

  Real evaluatePDF(const Real input) const {
    return ((input >= a_) ? ((input <= b_) ? coeff_*std::exp(-scale_*input) : 0.) : 0.);
  }

  Real evaluateCDF(const Real input) const {
    return ((input > a_) ? ((input < b_) ? (expa_-std::exp(-scale_*input))/diff_ : 1.) : 0.);
  }

  Real integrateCDF(const Real input) const {
    return ((input > a_) ? ((input < b_) ?
      (expa_*(input-a_) - (expa_ - std::exp(-scale_*input))/scale_)/diff_ :
        (expa_*(b_-a_) - (expa_ - expb_)/scale_)/diff_ + (input - b_)) : 0.);
  }

  Real invertCDF(const Real input) const {
    return ((input > 0.) ? ((input < 1.) ? -std::log(expa_-diff_*input)/scale_ : b_) : a_);
  }

  Real moment(const size_t m) const {
    Real val = 0., coeff = 0.;
    for (size_t i = 0; i < m+1; i++) {
      coeff = compute_coeff(m,i);
      val  += coeff*(std::pow(a_,i)*expa_-std::pow(b_,i)*expb_)/std::pow(scale_,m-i+1);
    } 
    return coeff_*val;
  }

  Real lowerBound(void) const {
    return a_;
  }
 
  Real upperBound(void) const {
    return b_;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 5;
    std::vector<Real> X(size,0.);
    std::vector<int> T(size,0);
    X[0] = a_-4.0*(Real)rand()/(Real)RAND_MAX;
    T[0] = 0;
    X[1] = a_; 
    T[1] = 1;
    X[2] = (b_-a_)*(Real)rand()/(Real)RAND_MAX + a_;
    T[2] = 0;
    X[3] = b_; 
    T[3] = 1;
    X[4] = b_+4.0*(Real)rand()/(Real)RAND_MAX;
    T[4] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
