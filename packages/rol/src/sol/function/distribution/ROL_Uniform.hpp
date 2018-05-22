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

#ifndef ROL_UNIFORM_HPP
#define ROL_UNIFORM_HPP

#include "ROL_Distribution.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class Uniform : public Distribution<Real> {
private:
  Real a_;
  Real b_;

public: 
  Uniform(const Real lo = 0., const Real up = 1.)
    : a_((lo < up) ? lo : up), b_((up > lo) ? up : lo) {}

  Uniform(ROL::ParameterList &parlist) {
    a_ = parlist.sublist("SOL").sublist("Distribution").sublist("Uniform").get("Lower Bound",0.);
    b_ = parlist.sublist("SOL").sublist("Distribution").sublist("Uniform").get("Upper Bound",1.);
    Real tmp = a_;
    a_ = std::min(a_,b_);
    b_ = std::max(b_,tmp);
  }

  Real evaluatePDF(const Real input) const {
    return ((input >= a_ && input <= b_) ? 1.0/(b_-a_) : 0.0);
  }

  Real evaluateCDF(const Real input) const {
    return ((input < a_) ? 0.0 : ((input > b_) ? 1.0 : (input-a_)/(b_-a_)));
  }

  Real integrateCDF(const Real input) const {
    return ((input < a_) ? 0.0 : ((input > b_) ? input - 0.5*(a_+b_) : 
              0.5*std::pow(input-a_,2.0)/(b_-a_)));
  }

  Real invertCDF(const Real input) const {
    return a_ + input*(b_-a_);
  }

  Real moment(const size_t m) const {
    return (std::pow(b_,m+1)-std::pow(a_,m+1))/((Real)(m+1)*(b_-a_));
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
    X[4] = b_+4.*(Real)rand()/(Real)RAND_MAX; 
    T[4] = 0;
    Distribution<Real>::test(X,T,outStream);
  }
};

}

#endif
