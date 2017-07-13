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

#ifndef ROL_GAMMA_HPP
#define ROL_GAMMA_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

#include <math.h>

namespace ROL {

template<class Real>
class Gamma : public Distribution<Real> {
private:
  Real shape_;
  Real scale_;
  Real gamma_shape_;
  Real coeff_;

  Real igamma(const Real s, const Real x) const {
    Real sum = 0., term = 1./s;
    size_t n = 1;
    while ( term != 0. ) {
      sum  += term;
      term *= x/(s+(Real)n);
      n++;
    }
    return std::pow(x,s)*std::exp(-x)*sum;
  }

public: 
  Gamma(const Real shape = 1., const Real scale = 1.)
    : shape_((shape > 0.) ? shape : 1.), scale_((scale > 0.) ? scale : 1.) {
    gamma_shape_ = tgamma(shape_);
    coeff_       = 1./(gamma_shape_*std::pow(scale_,shape_));
  }

  Gamma(Teuchos::ParameterList &parlist) {
    shape_ = parlist.sublist("SOL").sublist("Distribution").sublist("Gamma").get("Shape",1.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Gamma").get("Scale",1.);
    shape_ = (shape_ > 0.) ? shape_ : 1.;
    scale_ = (scale_ > 0.) ? scale_ : 1.;
    gamma_shape_ = tgamma(shape_);
    coeff_       = 1./(gamma_shape_*std::pow(scale_,shape_));
  }

  Real evaluatePDF(const Real input) const {
    return ((input <= 0.) ? 0. : coeff_*std::pow(input,shape_-1.)*std::exp(-input/scale_));
  }

  Real evaluateCDF(const Real input) const {
    return ((input <= 0.) ? 0. : igamma(shape_,input/scale_)/gamma_shape_);
  }

  Real integrateCDF(const Real input) const {
    Real x = input/scale_;
    return ((input <= 0.) ? 0. : (x*igamma(shape_,x) - igamma(shape_+1.,x))/scale_);
  }

  Real invertCDF(const Real input) const {
    if ( input <= 0. ) {
      return 0.;
    }
    Real x = input*gamma_shape_;
    Real fx = evaluateCDF(x)-input;
    Real s = 0., xs = 0., a = 1., tmp = 0.;
    for (size_t i = 0; i < 100; i++) {
      if ( std::abs(fx) < ROL_EPSILON<Real>() ) { break; }
      s   = -fx/evaluatePDF(x);
      a   = 1.0;
      xs  = x + a*s;
      tmp = fx;
      fx  = evaluateCDF(xs)-input;
      while ( std::abs(fx) > (1.0 - 1.e-4*a)*std::abs(tmp) ) {
        a *= 0.5;
        xs = x + a*s;
        fx = evaluateCDF(xs)-input;
      }
      x = xs;
    }
    return x;
  }

  Real moment(const size_t m) const {
    if ( m == 1 ) {
      return shape_*scale_;
    }
    if ( m == 2 ) {
      return shape_*scale_*scale_*(1. + shape_);
    }
    return std::pow(scale_,m)*tgamma(shape_+(Real)m)/gamma_shape_; 
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
    X[0] = -4.0*(Real)rand()/(Real)RAND_MAX;
    T[0] = 0;
    X[1] = 0.;
    T[1] = 1;
    X[2] = 4.0*(Real)rand()/(Real)RAND_MAX;
    T[2] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
