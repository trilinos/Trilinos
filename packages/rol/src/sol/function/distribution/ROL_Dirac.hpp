// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DIRAC_HPP
#define ROL_DIRAC_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Dirac : public Distribution<Real> {
private:
  Real data_;

public: 
  Dirac(const Real data = 0.) : data_(data) {}

  Dirac(ROL::ParameterList &parlist) {
    data_ = parlist.sublist("SOL").sublist("Distribution").sublist("Dirac").get("Location",0.);
  }

  Real evaluatePDF(const Real input) const {
    return ((input==data_) ? 1.0 : 0.0);
  }

  Real evaluateCDF(const Real input) const {
    return ((input >= data_) ? 1.0 : 0.0);
  }

  Real integrateCDF(const Real input) const {
    return ((input < data_) ? 0.0 : input);
  }

  Real invertCDF(const Real input) const {
    return data_;
  }

  Real moment(const size_t m) const {
    if (m==1) {
      return data_;
    }
    return std::pow(data_,(Real)m); 
  }

  Real lowerBound(void) const {
    return data_;
  }

  Real upperBound(void) const {
    return data_;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 0;
    std::vector<Real> X(size,4.*(Real)rand()/(Real)RAND_MAX - 2.);
    std::vector<int> T(size,0);
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
