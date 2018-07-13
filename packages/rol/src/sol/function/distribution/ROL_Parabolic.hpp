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

#ifndef ROL_PARABOLIC_HPP
#define ROL_PARABOLIC_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Parabolic : public Distribution<Real> {
private:
  Real a_;
  Real b_;

public: 
  Parabolic(const Real a = 0., const Real b = 1.)
    : a_(std::min(a,b)), b_(std::max(a,b)) {}

  Parabolic(ROL::ParameterList &parlist) {
    a_ = parlist.sublist("SOL").sublist("Distribution").sublist("Parabolic").get("Lower Bound",0.);
    b_ = parlist.sublist("SOL").sublist("Distribution").sublist("Parabolic").get("Upper Bound",1.);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);
  }

  Real evaluatePDF(const Real input) const {
    Real scale = 6.0/std::pow(b_-a_,3.0);
    return ((input >= a_ && input <= b_) ? scale*(input-a_)*(b_-input) : 0.);
  }

  Real evaluateCDF(const Real input) const {
    Real d1 = b_-a_, d2 = d1*d1, d3 = d2*d1;
    Real v1 = input-a_, v2 = v1*v1, v3 = v1*v2;
    return ((input < a_) ? 0. : ((input > b_) ? 1. : 
            3.0*v2/d2 - 2.0*v3/d3));
  }

  Real integrateCDF(const Real input) const {
    Real d0 = b_+a_, d1 = b_-a_, d2 = d1*d1, d3 = d2*d1;
    Real v1 = input-a_, v2 = v1*v1, v3 = v1*v2, v4 = v1*v3;
    return ((input < a_) ? 0. : 
           ((input > b_) ? input - 0.5*d0 :
             v3/d2 - 0.5*v4/d3));
  }

  Real invertCDF(const Real input) const {
    Real a  = a_-b_, b  = a_+b_, c  = 0.;
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
    Real p  = (Real)m;
    Real a1 = std::pow(a_,p+1), b1 = std::pow(b_,p+1);
    Real a2 = a1*a_,            b2 = b1*b_;
    Real a3 = a2*a_,            b3 = b2*b_;
    return 6./std::pow(b_-a_,3)
           * (-(b3-a3)/(p+3) + (a_+b_)*(b2-a2)/(p+2) - a_*b_*(b1-a1)/(p+1));
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
    X[0] = a_-4.0*(Real)rand()/(Real)RAND_MAX;
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
