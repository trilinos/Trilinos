// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LAPLACE_HPP
#define ROL_LAPLACE_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Laplace : public Distribution<Real> {
private:
  Real mean_;
  Real scale_;

  size_t compute_coeff(const size_t m, const size_t k) const {
    if ( k == 0 || m == 0 || m == 1 ) {
      return 1;
    }
    size_t val = 1;
    for (size_t i = m-k; i < m; i++) {
      val *= (i+1);
    }
    return val;
  }

public: 
  Laplace(const Real mean = 0., const Real scale = 1.)
    : mean_(mean), scale_(scale) {}

  Laplace(ROL::ParameterList &parlist) {
    mean_  = parlist.sublist("SOL").sublist("Distribution").sublist("Laplace").get("Mean",0.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Laplace").get("Scale",1.);
    scale_ = (scale_ > 0.) ? scale_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return 0.5*std::exp(-std::abs(input-mean_)/scale_)/scale_;
  }

  Real evaluateCDF(const Real input) const {
    return ((input < mean_) ? 0.5*std::exp((input-mean_)/scale_) : 
            1.-0.5*std::exp(-(input-mean_)/scale_));
  }

  Real integrateCDF(const Real input) const {
    return ((input < mean_) ? 0.5*scale_*std::exp((input-mean_)/scale_) : 
            (input-mean_)+0.5*scale_*std::exp(-(input-mean_)/scale_));
  }

  Real invertCDF(const Real input) const {
    Real sgn = ((input < 0.5) ? -1. : ((input > 0.5) ? 1. : 0.0));
    return mean_ - scale_*sgn*std::log(1.-2.*std::abs(input-0.5));
  }

  Real moment(const size_t m) const {
    if ( m == 1 ) {
      return mean_;
    }
    if ( m == 2 ) {
      return std::pow(mean_,2) + 2.*std::pow(scale_,2);
    }
    Real coeff = 0., val = 0.;
    for (size_t k = 0; k < m+1; k++) {
      if ( k%2 == 0 ) {
        coeff = compute_coeff(m,k);
        val  += coeff*std::pow(scale_,k)*std::pow(mean_,m-k);
      }
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
