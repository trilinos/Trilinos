// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BETA_HPP
#define ROL_BETA_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

#include <math.h>

namespace ROL {

template<class Real>
class Beta : public Distribution<Real> {
private:
  Real shape1_;
  Real shape2_;
  Real coeff_;

  std::vector<Real> pts_;
  std::vector<Real> wts_;

  void initializeQuadrature(void) {
    pts_.clear(); pts_.resize(20,0.); wts_.clear(); wts_.resize(20,0.);
    wts_[0]  = 0.1527533871307258; pts_[0]  = -0.0765265211334973;
    wts_[1]  = 0.1527533871307258; pts_[1]  =  0.0765265211334973;
    wts_[2]  = 0.1491729864726037; pts_[2]  = -0.2277858511416451;
    wts_[3]  = 0.1491729864726037; pts_[3]  =  0.2277858511416451;
    wts_[4]  = 0.1420961093183820; pts_[4]  = -0.3737060887154195;
    wts_[5]  = 0.1420961093183820; pts_[5]  =  0.3737060887154195;
    wts_[6]  = 0.1316886384491766; pts_[6]  = -0.5108670019508271;
    wts_[7]  = 0.1316886384491766; pts_[7]  =  0.5108670019508271;
    wts_[8]  = 0.1181945319615184; pts_[8]  = -0.6360536807265150;
    wts_[9]  = 0.1181945319615184; pts_[9]  =  0.6360536807265150;
    wts_[10] = 0.1019301198172404; pts_[10] = -0.7463319064601508;
    wts_[11] = 0.1019301198172404; pts_[11] =  0.7463319064601508;
    wts_[12] = 0.0832767415767048; pts_[12] = -0.8391169718222188;
    wts_[13] = 0.0832767415767048; pts_[13] =  0.8391169718222188;
    wts_[14] = 0.0626720483341091; pts_[14] = -0.9122344282513259;
    wts_[15] = 0.0626720483341091; pts_[15] =  0.9122344282513259;
    wts_[16] = 0.0406014298003869; pts_[16] = -0.9639719272779138;
    wts_[17] = 0.0406014298003869; pts_[17] =  0.9639719272779138;
    wts_[18] = 0.0176140071391521; pts_[18] = -0.9931285991850949;
    wts_[19] = 0.0176140071391521; pts_[19] =  0.9931285991850949;
    for (size_t i = 0; i < 20; i++) {
      wts_[i] *= 0.5;
      pts_[i] += 1.; pts_[i] *= 0.5;
    }
  }

  Real ibeta(const Real x) const {
    Real pt = 0., wt = 0., sum = 0.;
    for (size_t i = 0; i < pts_.size(); i++) {
      wt   = x*wts_[i];
      pt   = x*pts_[i];
      sum += wt*std::pow(pt,shape1_-1)*std::pow(1.-pt,shape2_-1);
    }
    return sum;
  }

public: 
  Beta(const Real shape1 = 2., const Real shape2 = 2.)
    : shape1_((shape1 > 0.) ? shape1 : 2.), shape2_((shape2 > 0.) ? shape2 : 2.) {
    coeff_ = tgamma(shape1_+shape2_)/(tgamma(shape1_)*tgamma(shape2_));
    initializeQuadrature();
  }

  Beta(ROL::ParameterList &parlist) {
    shape1_ = parlist.sublist("SOL").sublist("Distribution").sublist("Beta").get("Shape 1",2.);
    shape2_ = parlist.sublist("SOL").sublist("Distribution").sublist("Beta").get("Shape 2",2.);
    shape1_ = (shape1_ > 0.) ? shape1_ : 2.;
    shape2_ = (shape2_ > 0.) ? shape2_ : 2.;
    coeff_ = tgamma(shape1_+shape2_)/(tgamma(shape1_)*tgamma(shape2_));
    initializeQuadrature();
  }

  Real evaluatePDF(const Real input) const {
    return ((input > 0.) ? ((input < 1.) ?
             coeff_*std::pow(input,shape1_-1.)*std::pow(1.-input,shape2_-1) : 0.) : 0.);
  }

  Real evaluateCDF(const Real input) const {
    return ((input > 0.) ? ((input < 1.) ? coeff_*ibeta(input) : 1.) : 0.);
  }

  Real integrateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Beta): Beta integrateCDF not implemented!");
  }

  Real invertCDF(const Real input) const {
    if ( input <= ROL_EPSILON<Real>() ) {
      return 0.;
    }
    if ( input >= 1.-ROL_EPSILON<Real>() ) {
      return 1.;
    }
    Real a  = ROL_EPSILON<Real>(), b = 1.-ROL_EPSILON<Real>(), c  = 0.;
    Real fa = evaluateCDF(a) - input;
    Real fc = 0.;
    Real sa = ((fa < 0.) ? -1. : ((fa > 0.) ? 1. : 0.));
    Real sc = 0.;
    for (size_t i = 0; i < 100; i++) {
      c  = (a+b)*0.5;
      fc = evaluateCDF(c) - input;
      sc = ((fc < 0.) ? -1. : ((fc > 0.) ? 1. : 0.));
      if ( fc == 0. || (b-a)*0.5 < ROL_EPSILON<Real>() ) {
        break;
      }
      if ( sc == sa ) { a = c; fa = fc; sa = sc; }
      else            { b = c; }
    }
    return c;
  }

  Real moment(const size_t m) const {
    if ( m == 1 ) {
      return shape1_/(shape1_ + shape2_);
    }
    if ( m == 2 ) {
      return shape1_*(shape2_/(shape1_+shape2_+1.) + shape1_)/std::pow(shape1_+shape2_,2);
    }
    Real val = 1.;
    for (size_t i = 0; i < m; i++) {
      val *= (shape1_ + (Real)i)/(shape1_ + shape2_ + (Real)i);
    }
    return val;
  }

  Real lowerBound(void) const {
    return 0.;
  }

  Real upperBound(void) const {
    return 1.;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 5;
    std::vector<Real> X(size,0.);
    std::vector<int> T(size,0);
    X[0] = -4.*(Real)rand()/(Real)RAND_MAX;
    T[0] = 0;
    X[1] = 0.; 
    T[1] = 1;
    X[2] = 0.5*(Real)rand()/(Real)RAND_MAX;
    T[2] = 0;
    X[3] = 1.; 
    T[3] = 1;
    X[4] = 1.+4.*(Real)rand()/(Real)RAND_MAX;
    T[4] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
