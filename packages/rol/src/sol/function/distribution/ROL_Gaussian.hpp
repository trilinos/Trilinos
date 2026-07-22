// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_GAUSSIAN_HPP
#define ROL_GAUSSIAN_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Gaussian : public Distribution<Real> {
private:
  Real mean_;
  Real variance_;

  std::vector<Real> a_;
  std::vector<Real> b_;
  std::vector<Real> c_;
  std::vector<Real> d_;

  Real erfi(const Real p) const {
    const Real zero(0), half(0.5), one(1), two(2), pi(ROL::ScalarTraits<Real>::pi());
    Real val(0), z(0);
    if ( std::abs(p) > static_cast<Real>(0.7) ) {
      Real sgn = (p < zero) ? -one : one;
      z   = std::sqrt(-std::log((one-sgn*p)*half));
      val = sgn*(((c_[3]*z+c_[2])*z+c_[1])*z + c_[0])/((d_[1]*z+d_[0])*z + one);
    }
    else {
      z   = p*p;
      val = p*(((a_[3]*z+a_[2])*z+a_[1])*z + a_[0])/((((b_[3]*z+b_[2])*z+b_[1])*z+b_[0])*z+one);
    }
    val -= (erf(val)-p)/(two/std::sqrt(pi) * std::exp(-val*val));
    val -= (erf(val)-p)/(two/std::sqrt(pi) * std::exp(-val*val));
    return val;
  }

public: 

  Gaussian(const Real mean = 0., const Real variance = 1.)
    : mean_(mean), variance_((variance>0.) ? variance : 1.) {
    a_.clear(); a_.resize(4,0.); b_.clear(); b_.resize(4,0.);
    c_.clear(); c_.resize(4,0.); d_.clear(); d_.resize(2,0.);
    a_[0] =  0.886226899; a_[1] = -1.645349621; a_[2] =  0.914624893; a_[3] = -0.140543331;
    b_[0] = -2.118377725; b_[1] =  1.442710462; b_[2] = -0.329097515; b_[3] =  0.012229801;
    c_[0] = -1.970840454; c_[1] = -1.624906493; c_[2] =  3.429567803; c_[3] =  1.641345311;
    d_[0] =  3.543889200; d_[1] =  1.637067800;
  }

  Gaussian(ROL::ParameterList &parlist) {
    mean_     = parlist.sublist("SOL").sublist("Distribution").sublist("Gaussian").get("Mean",0.);
    variance_ = parlist.sublist("SOL").sublist("Distribution").sublist("Gaussian").get("Variance",1.);
    variance_ = (variance_ > 0.) ? variance_ : 1.;
    a_.clear(); a_.resize(4,0.); b_.clear(); b_.resize(4,0.);
    c_.clear(); c_.resize(4,0.); d_.clear(); d_.resize(2,0.);
    a_[0] =  0.886226899; a_[1] = -1.645349621; a_[2] =  0.914624893; a_[3] = -0.140543331;
    b_[0] = -2.118377725; b_[1] =  1.442710462; b_[2] = -0.329097515; b_[3] =  0.012229801;
    c_[0] = -1.970840454; c_[1] = -1.624906493; c_[2] =  3.429567803; c_[3] =  1.641345311;
    d_[0] =  3.543889200; d_[1] =  1.637067800;
  }

  Real evaluatePDF(const Real input) const {
    return std::exp(-std::pow(input-mean_,2)/(2.*variance_))/(std::sqrt(2.*ROL::ScalarTraits<Real>::pi()*variance_));
  }

  Real evaluateCDF(const Real input) const {
    const Real half(0.5), one(1), two(2);
    return half*(one+erf((input-mean_)/std::sqrt(two*variance_)));
  }

  Real integrateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Gaussian): Gaussian integrateCDF not implemented!");
  }

  Real invertCDF(const Real input) const {
    //return std::sqrt(2.*variance_)*erfi(2.*input-1.) + mean_;
    const Real zero(0), half(0.5), one(1), eight(8);
    const Real dev(std::sqrt(variance_)), eps(1.24419211485e-15);
    // Set lower and upper bounds to the mean plus/minus 8 standard
    //   -- deviations this ensures that 1-eps probability mass is
    //   -- covered by the interval.
    const Real lVal = mean_ - eight*dev;
    const Real uVal = mean_ + eight*dev;
    // If the input is outside of the interval (half*eps,1-half*eps)
    //   -- then set the return value to be either the lower or
    //   -- upper bound.  This case can occur with probability eps.
    if ( input <= half*eps )     { return lVal; }
    if ( input >= one-half*eps ) { return uVal; }
    // Determine maximum number of iterations.
    //   -- maxit is set to the number of iterations required to
    //   -- ensure that |b-a| < eps after maxit iterations.
    size_t maxit = static_cast<size_t>(one-std::log2(eps/(eight*dev)));
    maxit = (maxit < 1 ? 100 : maxit);
    // Run bisection to solve CDF(x) = input.
    Real a  = (input < half ? lVal  : mean_);
    Real b  = (input < half ? mean_ : uVal );
    Real c  = half*(a+b);
    Real fa = evaluateCDF(a) - input;
    Real fc = evaluateCDF(c) - input;
    Real sa = ((fa < zero) ? -one : ((fa > zero) ? one : zero));
    Real sc = ((fc < zero) ? -one : ((fc > zero) ? one : zero));
    for (size_t i = 0; i < maxit; ++i) {
      if ( std::abs(fc) < eps || (b-a)*half < eps ) {
        break;
      }
      if ( sc == sa ) { a = c; fa = fc; sa = sc; }
      else            { b = c; }
      // Compute interval midpoint.
      c  = (a+b)*half;
      fc = evaluateCDF(c) - input;
      sc = ((fc < zero) ? -one : ((fc > zero) ? one : zero));
    }
    return c;
  }

  Real moment(const size_t m) const {
    Real val = 0.;
    switch(m) {
      case 1: val = mean_;                                         break;
      case 2: val = std::pow(mean_,2) + variance_;                 break;
      case 3: val = std::pow(mean_,3)
                    + 3.*mean_*variance_;                          break;
      case 4: val = std::pow(mean_,4)
                    + 6.*std::pow(mean_,2)*variance_
                    + 3.*std::pow(variance_,2);                    break;
      case 5: val = std::pow(mean_,5)
                    + 10.*std::pow(mean_,3)*variance_
                    + 15.*mean_*std::pow(variance_,2);             break;
      case 6: val = std::pow(mean_,6)
                    + 15.*std::pow(mean_,4)*variance_
                    + 45.*std::pow(mean_*variance_,2)
                    + 15.*std::pow(variance_,3);                   break;
      case 7: val = std::pow(mean_,7)
                    + 21.*std::pow(mean_,5)*variance_
                    + 105.*std::pow(mean_,3)*std::pow(variance_,2)
                    + 105.*mean_*std::pow(variance_,3);            break;
      case 8: val = std::pow(mean_,8)
                    + 28.*std::pow(mean_,6)*variance_
                    + 210.*std::pow(mean_,4)*std::pow(variance_,2)
                    + 420.*std::pow(mean_,2)*std::pow(variance_,3)
                    + 105.*std::pow(variance_,4);                  break;
      default:
        ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
          ">>> ERROR (ROL::Distribution): Gaussian moment not implemented for m > 8!");
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
