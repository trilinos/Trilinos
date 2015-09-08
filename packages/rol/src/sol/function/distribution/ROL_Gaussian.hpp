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

#ifndef ROL_GAUSSIAN_HPP
#define ROL_GAUSSIAN_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Gaussian : public Distribution<Real> {
private:
  Real mean_;
  Real variance_;

public: 

  Gaussian(const Real mean = 0., const Real variance = 1.)
    : mean_(mean), variance_((variance>0.) ? variance : 1.) {}

  Gaussian(Teuchos::ParameterList &parlist) {
    mean_     = parlist.sublist("SOL").sublist("Distribution").sublist("Gaussian").get("Mean",0.);
    variance_ = parlist.sublist("SOL").sublist("Distribution").sublist("Gaussian").get("Variance",1.);
    variance_ = (variance_ > 0.) ? variance_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return std::exp(-std::pow(input-mean_,2.0)/(2.0*variance_))/(std::sqrt(2.0*M_PI*variance_));
  }

  Real evaluateCDF(const Real input) const {
    return 0.5*(1.0+erf((input-mean_)/std::sqrt(2.0*variance_)));
  }

  Real integrateCDF(const Real input) const {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Gaussian): Gaussian integrateCDF not implemented!");
    return ((input < mean_) ? 0.0 : input);
  }

  Real invertCDF(const Real input) const {
    std::vector<Real> coeff;
    Real x   = 2.0*input - 1.0;
    Real c   = 1.0;
    Real tmp = c * (std::sqrt(M_PI)/2.0 * x);
    Real val = tmp;
    coeff.push_back(c);
    int  cnt = 1;
    while (std::abs(tmp) > std::sqrt(ROL_EPSILON)*std::abs(val)) {
      c = 0.0;
      for ( unsigned i = 0; i < coeff.size(); i++ ) {
        c += coeff[i]*coeff[coeff.size()-1-i]/((i+1)*(2*i+1));
      }
      tmp  = c/(2.0*(Real)cnt+1.0) * std::pow(std::sqrt(M_PI)/2.0 * x,2.0*(Real)cnt+1.0);
      val += tmp;
      coeff.push_back(c);
      cnt++;
    }
    return std::sqrt(2*variance_)*val + mean_;
  }

  Real moment(const size_t m) const {
    Real val = 0.0;
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
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
          ">>> ERROR (ROL::Distribution): Gaussian moment not implemented for m > 8!");
    }
    return val;
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
