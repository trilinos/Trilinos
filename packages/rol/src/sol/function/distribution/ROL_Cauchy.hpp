// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CAUCHY_HPP
#define ROL_CAUCHY_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Cauchy : public Distribution<Real> {
private:
  Real loc_;
  Real scale_;

public: 
  Cauchy(const Real loc = 0., const Real scale = 1.)
    : loc_(loc), scale_((scale>0.) ? scale : 1.) {}

  Cauchy(ROL::ParameterList &parlist) {
    loc_   = parlist.sublist("SOL").sublist("Distribution").sublist("Cauchy").get("Location",0.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Cauchy").get("Scale",1.);
    scale_ = (scale_>0.) ? scale_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return 1./(ROL::ScalarTraits<Real>::pi()*scale_*(1.+std::pow((input-loc_)/scale_,2.)));
  }

  Real evaluateCDF(const Real input) const {
    return 0.5+atan((input-loc_)/scale_)/ROL::ScalarTraits<Real>::pi();
  }

  Real integrateCDF(const Real input) const {
    Real v = input-loc_;
    return 0.5*input + (v*atan(v/scale_) - 0.5*scale_*std::log(v*v+scale_*scale_))/ROL::ScalarTraits<Real>::pi();
  }

  Real invertCDF(const Real input) const {
    return loc_+scale_*tan(ROL::ScalarTraits<Real>::pi()*(input-0.5));
  }

  Real moment(const size_t m) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Cauchy): Cauchy moments are undefined!");
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
