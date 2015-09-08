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
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Exponential : public Distribution<Real> {
private:
  Real loc_;
  Real scale_;

  size_t factorial(const size_t m) const {
    return (m==1 ? m : m * factorial(m-1));
  }

public: 
  Exponential(const Real loc = 0., const Real scale = 1.)
    : loc_(loc), scale_((scale>0.) ? scale : 1.) {}

  Exponential(Teuchos::ParameterList &parlist) {
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
    Real val = 0., coeff = 0., factm = (Real)factorial(m);
    for (size_t i = 0; i < m; i++) {
      coeff = factm/(Real)factorial(m-i);
      val  += coeff*std::pow(loc_,(Real)(m-i))/std::pow(scale_,(Real)i);
    } 
    return val;
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
