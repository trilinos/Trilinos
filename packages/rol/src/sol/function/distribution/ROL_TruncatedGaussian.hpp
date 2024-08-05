// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUNCATEDGAUSSIAN_HPP
#define ROL_TRUNCATEDGAUSSIAN_HPP

#include "ROL_Distribution.hpp"
#include "ROL_Gaussian.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class TruncatedGaussian : public Distribution<Real> {
private:
  Real a_;
  Real b_;
  Real mean_;
  Real sdev_;

  ROL::Ptr<Gaussian<Real> > gauss_;

  Real alpha_;
  Real beta_;
  Real phi_;
  Real Z_;

public: 

  TruncatedGaussian(const Real lo = -1., const Real up = 1.,
                    const Real mean = 0., const Real sdev = 1.)
    : a_((lo < up) ? lo : up), b_((up > lo) ? up : lo),
      mean_(mean), sdev_((sdev>0) ? sdev : 1) {
    //Real var = sdev_*sdev_;
    gauss_ = ROL::makePtr<Gaussian<Real>>();
    alpha_ = (a_-mean_)/sdev_;
    beta_  = (b_-mean_)/sdev_;
    phi_   = gauss_->evaluateCDF(alpha_);
    Z_     = gauss_->evaluateCDF(beta_)-gauss_->evaluateCDF(alpha_);
  }

  TruncatedGaussian(ROL::ParameterList &parlist) {
    const Real zero(0), one(1);
    ROL::ParameterList TGlist
      = parlist.sublist("SOL").sublist("Distribution").sublist("Truncated Gaussian");
    a_ = TGlist.get("Lower Bound",-one);
    b_ = TGlist.get("Upper Bound", one);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);

    mean_ = TGlist.get("Mean",zero);
    sdev_ = TGlist.get("Standard Deviation",one);
    sdev_ = (sdev_ > zero) ? sdev_ : one;

    //Real var = sdev_*sdev_;
    gauss_ = ROL::makePtr<Gaussian<Real>>();
    alpha_ = (a_-mean_)/sdev_;
    beta_  = (b_-mean_)/sdev_;
    phi_   = gauss_->evaluateCDF(alpha_);
    Z_     = gauss_->evaluateCDF(beta_)-gauss_->evaluateCDF(alpha_);
  }

  Real evaluatePDF(const Real input) const {
    const Real zero(0), xi = (input-mean_)/sdev_;
    return ((input <= a_) ? zero : ((input >= b_) ? zero :
             gauss_->evaluatePDF(xi)/(sdev_*Z_)));
  }

  Real evaluateCDF(const Real input) const {
    const Real zero(0), one(1), xi = (input-mean_)/sdev_;
    return ((input <= a_) ? zero : ((input >= b_) ? one : 
             (gauss_->evaluateCDF(xi)-phi_)/Z_));
  }

  Real integrateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::TruncatedGaussian): Truncated Gaussian integrateCDF not implemented!");
    //return ((input < 0.5*(a_+b_)) ? 0.0 : input - 0.5*(a_+b_));
  }

  Real invertCDF(const Real input) const {
    const Real x = gauss_->invertCDF(Z_*input+phi_);
    return sdev_*x + mean_;
  }

  Real moment(const size_t m) const {
    const Real phiA  = gauss_->evaluatePDF(alpha_);
    const Real phiB  = gauss_->evaluatePDF(beta_);
    const Real mean  = mean_ + sdev_*(phiA-phiB)/Z_;
    const Real var   = sdev_*sdev_;
    const Real one(1);
    Real val(0);
    switch(m) {
      case 1: val = mean;                                                                       break;
      case 2: val = var*(one+(alpha_*phiA-beta_*phiB)/Z_-std::pow((phiA-phiB)/Z_,2))+mean*mean; break;
      default:
        ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
          ">>> ERROR (ROL::TruncatedGaussian): Truncated Gaussian moment not implemented for m > 2!");
    }
    return val;
  }

  Real lowerBound(void) const {
    return a_;
  }
 
  Real upperBound(void) const {
    return b_;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 5;
    std::vector<Real> X(size,0);
    std::vector<int> T(size,0);
    const Real four(4);
    X[0] = a_-four*(Real)rand()/(Real)RAND_MAX; 
    T[0] = 0;
    X[1] = a_; 
    T[1] = 1;
    X[2] = (b_-a_)*(Real)rand()/(Real)RAND_MAX + a_; 
    T[2] = 0;
    X[3] = b_; 
    T[3] = 1;
    X[4] = b_+four*(Real)rand()/(Real)RAND_MAX; 
    T[4] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
