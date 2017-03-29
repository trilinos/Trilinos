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

#ifndef ROL_RAISEDCOSINE_HPP
#define ROL_RAISEDCOSINE_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class RaisedCosine : public Distribution<Real> {
private:
  Real mean_;
  Real var_;

  size_t factorial(const size_t m) const {
    return (m==1 ? m : m * factorial(m-1));
  }

public: 
  RaisedCosine(const Real mean = 0.5, const Real var = 0.5)
    : mean_(mean), var_(((var>0.) ? var : 0.5)) {}

  RaisedCosine(Teuchos::ParameterList &parlist) {
    mean_ = parlist.sublist("SOL").sublist("Distribution").sublist("Raised Cosine").get("Mean",0.5);
    var_  = parlist.sublist("SOL").sublist("Distribution").sublist("Raised Cosine").get("Scale",0.5);
    var_  = (var_ > 0.) ? var_ : 0.5;
  }

  Real evaluatePDF(const Real input) const {
    Real a = mean_-var_, b = mean_+var_;
    return ((input >= a && input <= b) ?
             (1.+std::cos(Teuchos::ScalarTraits<Real>::pi()*(input-mean_)/var_))/(2.0*var_) : 0.);
  }

  Real evaluateCDF(const Real input) const {
    Real a = mean_-var_, b = mean_+var_;
    return ((input < a) ? 0. : ((input > b) ? 1. : 
            0.5*(1.+(input-mean_)/var_+std::sin(Teuchos::ScalarTraits<Real>::pi()*(input-mean_)/var_)/Teuchos::ScalarTraits<Real>::pi())));
  }
  Real integrateCDF(const Real input) const {
    Real a = mean_-var_, b = mean_+var_;
    Real v = input-mean_;
    return ((input < a) ? 0. : ((input > b) ? input - var_ : 
            0.5*(v+0.5*v*v/var_-var_*((std::cos(Teuchos::ScalarTraits<Real>::pi()*v/var_)+1.) /
                                      (Teuchos::ScalarTraits<Real>::pi()*Teuchos::ScalarTraits<Real>::pi())-0.5))));
  }
  Real invertCDF(const Real input) const {
    Real a = mean_-var_, b = mean_+var_, c  = 0.;
    Real fa = evaluateCDF(a) - input;
    Real fc = 0.;
    Real sa = ((fa < 0.) ? -1. : ((fa > 0.) ? 1. : 0.));
    Real sc = 0.;
    for (size_t i = 0; i < 100; i++) {
      c  = (a+b)*0.5;
      fc = evaluateCDF(c) - input;
      sc = ((fc < 0.) ? -1. : ((fc > 0.) ? 1. : 0.));
      if ( fc == 0. || (b-a)*0.5 < ROL_EPSILON<Real>() ) {
        break;
      }
      if ( sc == sa ) { a = c; fa = fc; sa = sc; }
      else            { b = c; }
    } 
    return c;
  }

  Real moment(const size_t m) const {
    Real a = mean_-var_, b = mean_+var_;
    Real am = std::pow(a,m+1), bm = std::pow(b,m+1);
    Real omega = Teuchos::ScalarTraits<Real>::pi()/var_, phi = -Teuchos::ScalarTraits<Real>::pi()*mean_/var_;
    Real val_cos = 0., val_sin = 0.;
    for (size_t k = 0; k < (m-1)/2; k++) {
      val_cos += ((k%2==0) ? 1. : -1.)*factorial(m)/(factorial(m-2*k-1)*std::pow(omega,2+2*k))
                *(std::pow(b,m-2*k-1)*std::cos(omega*b+phi)-std::pow(a,m-2*k-1)*std::cos(omega*a+phi));
    }
    for (size_t k = 0; k < m/2; k++) {
      val_sin += ((k%2==0) ? 1. : -1.)*factorial(m)/(factorial(m-2*k)*std::pow(omega,1+2*k))
                *(std::pow(b,m-2*k)*std::sin(omega*b+phi)-std::pow(a,m-2*k)*std::sin(omega*a+phi));
    }
    return 0.5*((bm-am)/((Real)m+1) + val_cos + val_sin)/var_;
  }

  Real lowerBound(void) const {
    return mean_-var_;
  }
 
  Real upperBound(void) const {
    return mean_+var_;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 5;
    std::vector<Real> X(size,0.);
    std::vector<int> T(size,0);
    X[0] = mean_-var_-4.*(Real)rand()/(Real)RAND_MAX; 
    T[0] = 0;
    X[1] = mean_-var_; 
    T[1] = 1;
    X[2] = (2.*var_)*(Real)rand()/(Real)RAND_MAX + (mean_-var_); 
    T[2] = 0;
    X[3] = mean_+var_; 
    T[3] = 1;
    X[4] = mean_+var_+4.*(Real)rand()/(Real)RAND_MAX; 
    T[4] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
