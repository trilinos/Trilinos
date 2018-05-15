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

#ifndef ROL_KUMARASWAMY_HPP
#define ROL_KUMARASWAMY_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

#include <math.h>

namespace ROL {

template<class Real>
class Kumaraswamy : public Distribution<Real> {
private:
  Real a_;
  Real b_;
  Real exp1_;
  Real exp2_;
  Real bGAMMAb_;

  size_t nchoosek(const size_t n, const size_t k) const {
    return ((k==0) ? 1 : (n*nchoosek(n-1,k-1))/k);
  }

public: 
  Kumaraswamy(const Real a = 0., const Real b = 1.,
              const Real exp1 = 0.5, const Real exp2 = 0.5)
    : a_(std::min(a,b)), b_(std::max(a,b)),
      exp1_((exp1>0.) ? exp1 : 0.5), exp2_((exp2>0.) ? exp2 : 0.5) {
    bGAMMAb_ = exp2_*tgamma(exp2_);
  }

  Kumaraswamy(ROL::ParameterList &parlist) {
    a_ = parlist.sublist("SOL").sublist("Distribution").sublist("Kumaraswamy").get("Lower Bound",0.);
    b_ = parlist.sublist("SOL").sublist("Distribution").sublist("Kumaraswamy").get("Upper Bound",1.);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);

    exp1_ = parlist.sublist("SOL").sublist("Distribution").sublist("Kumaraswamy").get("Exponent 1",0.5);
    exp2_ = parlist.sublist("SOL").sublist("Distribution").sublist("Kumaraswamy").get("Exponent 2",0.5);
    exp1_ = (exp1_ > 0.) ? exp1_ : 0.5;
    exp2_ = (exp2_ > 0.) ? exp2_ : 0.5;

    bGAMMAb_ = exp2_*tgamma(exp2_);
  }

  Real evaluatePDF(const Real input) const {
    Real x = (input - a_)/(b_-a_);
    return ((x <= 0.) ? 0. : ((x >= 1.) ? 0. : 
             exp1_*exp2_*std::pow(x,exp1_-1)*std::pow(1.-std::pow(x,exp1_),exp2_-1)));
  }

  Real evaluateCDF(const Real input) const {
    Real x = (input - a_)/(b_-a_);
    return ((x <= 0.) ? 0. : ((x >= 1.) ? 1. : 1.-std::pow(1.-std::pow(x,exp1_),exp2_)));
  }

  Real integrateCDF(const Real input) const {
    ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Kumaraswamy): Kumaraswamy integrateCDF not implemented!");
    return ((input < 0.5*(a_+b_)) ? 0.0 : input - 0.5*(a_+b_));
  }

  Real invertCDF(const Real input) const {
    Real x = std::pow(1.-std::pow(1.-input,1./exp2_),1./exp1_);
    return x*(b_-a_) + a_;
  }

  Real moment(const size_t m) const {
    Real val = 0., binom = 0., moment = 0.;
    for (size_t i = 0; i < m+1; i++) {
      moment = bGAMMAb_*tgamma(1.+(Real)i/exp1_)/tgamma(1.+exp2_+(Real)i/exp1_);
      binom  = (Real)nchoosek(m,i);
      val   += binom*std::pow(a_,m-i)*std::pow(b_-a_,i+1)*moment;
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
    std::vector<Real> X(size,0.);
    std::vector<int> T(size,0);
    X[0] = a_-4.*(Real)rand()/(Real)RAND_MAX; 
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
