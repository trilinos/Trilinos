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

#ifndef ROL_CAUCHY_HPP
#define ROL_CAUCHY_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Cauchy : public Distribution<Real> {
private:
  Real loc_;
  Real scale_;

public: 
  Cauchy(const Real loc = 0., const Real scale = 1.)
    : loc_(loc), scale_((scale>0.) ? scale : 1.) {}

  Cauchy(Teuchos::ParameterList &parlist) {
    loc_   = parlist.sublist("SOL").sublist("Distribution").sublist("Cauchy").get("Location",0.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Cauchy").get("Scale",1.);
    scale_ = (scale_>0.) ? scale_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return 1./(Teuchos::ScalarTraits<Real>::pi()*scale_*(1.+std::pow((input-loc_)/scale_,2.)));
  }

  Real evaluateCDF(const Real input) const {
    return 0.5+atan((input-loc_)/scale_)/Teuchos::ScalarTraits<Real>::pi();
  }

  Real integrateCDF(const Real input) const {
    Real v = input-loc_;
    return 0.5*input + (v*atan(v/scale_) - 0.5*scale_*std::log(v*v+scale_*scale_))/Teuchos::ScalarTraits<Real>::pi();
  }

  Real invertCDF(const Real input) const {
    return loc_+scale_*tan(Teuchos::ScalarTraits<Real>::pi()*(input-0.5));
  }

  Real moment(const size_t m) const {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
      ">>> ERROR (ROL::Cauchy): Cauchy moments are undefined!");
    return 0.;
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
