// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LOGISTIC_HPP
#define ROL_LOGISTIC_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Logistic : public Distribution<Real> {
private:
  Real mean_;
  Real var_;

public: 
  Logistic(const Real mean = 0., const Real var = 1.)
    : mean_(mean), var_((var>0.) ? var : 1.) {}

  Logistic(ROL::ParameterList &parlist) {
    mean_ = parlist.sublist("SOL").sublist("Distribution").sublist("Logistic").get("Mean",0.);
    var_  = parlist.sublist("SOL").sublist("Distribution").sublist("Logistic").get("Scale",1.);
    var_  = (var_ > 0.) ? var_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    Real val = std::exp(-(input-mean_)/var_);
    return val/(var_*std::pow(1.0+val,2.0));
  }

  Real evaluateCDF(const Real input) const {
    Real val = std::exp(-(input-mean_)/var_);
    return 1.0/(1.0+val);
  }

  Real integrateCDF(const Real input) const {
    Real val = std::exp(-(input-mean_)/var_);
    return (input-mean_) + var_*std::log(1.0+val);
  }

  Real invertCDF(const Real input) const {
    return mean_ + var_*std::log(input/(1.0-input));
  }

  Real moment(const size_t m) const {
    Real val = 0.;
    switch(m) {
      case 1: val = mean_;                                        break;
      case 2: val = std::pow(var_*ROL::ScalarTraits<Real>::pi(),2)/3. + std::pow(mean_,2); break;
      default:
        ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
          ">>> ERROR (ROL::Logistic): Logistic moment not implemented for m > 2!");
    }
    return val;
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
