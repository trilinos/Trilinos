// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SMALE_HPP
#define ROL_SMALE_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Smale : public Distribution<Real> {
private:
  Real a_;
  Real b_;

public: 
  Smale(const Real a = 0., const Real b = 1.)
    : a_(std::min(a,b)), b_(std::max(a,b)) {}

  Smale(ROL::ParameterList &parlist) {
    a_ = parlist.sublist("SOL").sublist("Distribution").sublist("Smale").get("Lower Bound",0.);
    b_ = parlist.sublist("SOL").sublist("Distribution").sublist("Smale").get("Upper Bound",1.);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);
  }

  Real evaluatePDF(const Real input) const {
    Real val  = std::pow(input-a_,2)+4.*b_*b_;
    Real root = std::sqrt(val);
    return 2.0*b_*b_/(val*root);
  }

  Real evaluateCDF(const Real input) const {
    Real val  = std::pow(input-a_,2)+4.*b_*b_;
    Real root = std::sqrt(val);
    return 0.5*(1.0+input/root);
  }

  Real integrateCDF(const Real input) const {
    Real val  = std::pow(input-a_,2)+4.*b_*b_;
    Real root = std::sqrt(val);
    return 0.5*(input+root);
  }

  Real invertCDF(const Real input) const {
    Real x   = a_;
    Real fx  = evaluateCDF(x)-input;
    Real s   = 0.0;
    Real xs  = 0.0;
    Real a   = 1.0;
    Real tmp = 0.0;
    for (int i = 0; i < 100; i++) {
      if ( std::abs(fx) < ROL_EPSILON<Real>() ) { break; }
      s   = -fx/evaluatePDF(x);
      a   = 1.0;
      xs  = x + a*s;
      tmp = fx;
      fx  = evaluateCDF(xs)-input;
      while ( std::abs(fx) > (1.0 - 1.e-4*a)*std::abs(tmp) ) {
        a *= 0.5;
        xs = x + a*s;
        fx = evaluateCDF(xs)-input;
      }
      x = xs;
    }
    return x;
  }

  Real moment(const size_t m) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Smale): Smale moment is not implemented!");
  }

  Real lowerBound(void) const {
    return ROL_NINF<Real>();
  }
 
  Real upperBound(void) const {
    return ROL_INF<Real>();
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 1;
    std::vector<Real> X(size,4.*(Real)rand()/(Real)RAND_MAX - 2.);
    std::vector<int> T(size,0);
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
