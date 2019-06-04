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

#ifndef ROL_EXPONENTIAL_HPP
#define ROL_EXPONENTIAL_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Exponential : public Distribution<Real> {
private:
  Real loc_;
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
  Exponential(const Real loc = 0., const Real scale = 1.)
    : loc_(loc), scale_((scale>0.) ? scale : 1.) {}

  Exponential(ROL::ParameterList &parlist) {
    loc_   = parlist.sublist("SOL").sublist("Distribution").sublist("Exponential").get("Location",0.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Exponential").get("Scale",1.);
    scale_ = (scale_ > 0.) ? scale_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return ((input >= loc_) ? scale_*std::exp(-scale_*(input-loc_)) : 0.);
  }

  Real evaluateCDF(const Real input) const {
    return ((input >= loc_) ? 1.-std::exp(-scale_*(input-loc_)) : 0.);
  }

  Real integrateCDF(const Real input) const {
    return ((input >= loc_) ?
      (input-loc_) - (1.-std::exp(-scale_*(input-loc_)))/scale_ : 0.);
  }

  Real invertCDF(const Real input) const {
    return -std::log(1.-input)/scale_;
  }

  Real moment(const size_t m) const {
    Real val = 0., coeff = 0.;
    for (size_t i = 0; i < m+1; i++) {
      coeff = compute_coeff(m,i);
      val  += coeff*std::pow(loc_,(Real)(m-i))/std::pow(scale_,(Real)i);
    } 
    return val;
  }

  Real lowerBound(void) const {
    return 0.;
  }

  Real upperBound(void) const {
    return ROL_INF<Real>();
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 3;
    std::vector<Real> X(size,0.);
    std::vector<int> T(size,0);
    X[0] = loc_-4.0*(Real)rand()/(Real)RAND_MAX;
    T[0] = 0;
    X[1] = loc_;
    T[1] = 1;
    X[2] = loc_+4.0*(Real)rand()/(Real)RAND_MAX;
    T[2] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
