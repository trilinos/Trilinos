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

#ifndef ROL_TRIANGLE_HPP
#define ROL_TRIANGLE_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Triangle : public Distribution<Real> {
private:
  Real a_;
  Real b_;
  Real c_;

public: 
  Triangle(const Real a = 0., const Real b = 0.5, const Real c = 1.)
    : a_(std::min(a,std::min(b,c))),
      b_(std::max(std::min(a,b),std::min(std::max(a,b),c))),
      c_(std::max(a,std::max(b,c))) {}

  Triangle(Teuchos::ParameterList &parlist) {
    Real a = parlist.sublist("SOL").sublist("Distribution").sublist("Triangle").get("Lower Bound",0.);
    Real b = parlist.sublist("SOL").sublist("Distribution").sublist("Triangle").get("Peak Location",0.5);
    Real c = parlist.sublist("SOL").sublist("Distribution").sublist("Triangle").get("Upper Bound",1.);
    a_ = std::min(a,std::min(b,c));
    b_ = std::max(std::min(a,b),std::min(std::max(a,b),c));
    c_ = std::max(a,std::max(b,c));
  }

  Real evaluatePDF(const Real input) const {
    Real d1 = b_-a_, d2 = c_-b_, d = c_-a_;
    return ((input >= a_ && input < b_) ? 2.0*(input-a_)/(d*d1) :
           ((input >= b_ && input < c_) ? 2.0*(c_-input)/(d*d2) : 
             0.0));
  }

  Real evaluateCDF(const Real input) const {
    Real d1 = b_-a_, d2 = c_-b_, d = c_-a_;
    return ((input < a_) ? 0.0 : 
           ((input >= a_ && input < b_) ? 
             std::pow(input-a_,2.0)/(d*d1) :
           ((input >= b_ && input < c_) ? 
             1.0-std::pow(c_-input,2.0)/(d*d2) : 
             1.0)));
  }

  Real integrateCDF(const Real input) const {
    Real d1 = b_-a_, d2 = c_-b_, d = c_-a_;
    return ((input < a_) ? 0.0 : 
           ((input >= a_ && input < b_) ? 
             std::pow(input-a_,3.0)/(3.0*d*d1) : 
           ((input >= b_ && input < c_) ?
             d1*d1/(3.0*d)+(input-b_)+(std::pow(c_-input,3.0)-d2*d2*d2)/(3.0*d*d2) :
             d1*d1/(3.0*d)+(input-b_)-d2*d2/(3.0*d))));
  }

  Real invertCDF(const Real input) const {
    Real d1 = b_-a_, d2 = c_-b_, d = c_-a_;
    return ((input <= d1/d) ? a_ + std::sqrt(input*d1*d) :
             c_ - std::sqrt((1.0-input)*d2*d));
  }

  Real moment(const size_t m) const {
    Real d1 = b_-a_, d2 = c_-b_, d = c_-a_;
    Real am1 = std::pow(a_,m+1), am2 = a_*am1;
    Real bm1 = std::pow(b_,m+1), bm2 = b_*bm1;
    Real cm1 = std::pow(c_,m+1), cm2 = c_*cm1;
    return (2./d)*(((bm2-am2)/((Real)m+2)-a_*(bm1-am1)/((Real)m+1))/d1
                  +(c_*(cm1-bm1)/((Real)m+1)-(cm2-bm2)/((Real)m+2))/d2);
  }

  Real lowerBound(void) const {
    return a_;
  }
 
  Real upperBound(void) const {
    return c_;
  }
 
  void test(std::ostream &outStream = std::cout ) const {
    size_t size = 7;
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
    X[4] = (c_-b_)*(Real)rand()/(Real)RAND_MAX + b_;
    T[4] = 0;
    X[5] = c_;
    T[5] = 1;
    X[6] = c_+4.*(Real)rand()/(Real)RAND_MAX;
    T[6] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
