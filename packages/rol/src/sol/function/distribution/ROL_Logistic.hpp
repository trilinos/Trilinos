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

#ifndef ROL_LOGISTIC_HPP
#define ROL_LOGISTIC_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Logistic : public Distribution<Real> {
private:
  Real mean_;
  Real var_;

public: 
  Logistic(const Real mean = 0., const Real var = 1.)
    : mean_(mean), var_((var>0.) ? var : 1.) {}

  Logistic(Teuchos::ParameterList &parlist) {
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
      case 2: val = std::pow(var_*Teuchos::ScalarTraits<Real>::pi(),2)/3. + std::pow(mean_,2); break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
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
