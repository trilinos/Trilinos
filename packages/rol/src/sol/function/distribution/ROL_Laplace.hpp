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

#ifndef ROL_LAPLACE_HPP
#define ROL_LAPLACE_HPP

#include "ROL_Distribution.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class Laplace : public Distribution<Real> {
private:
  Real mean_;
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
  Laplace(const Real mean = 0., const Real scale = 1.)
    : mean_(mean), scale_(scale) {}

  Laplace(Teuchos::ParameterList &parlist) {
    mean_  = parlist.sublist("SOL").sublist("Distribution").sublist("Laplace").get("Mean",0.);
    scale_ = parlist.sublist("SOL").sublist("Distribution").sublist("Laplace").get("Scale",1.);
    scale_ = (scale_ > 0.) ? scale_ : 1.;
  }

  Real evaluatePDF(const Real input) const {
    return 0.5*std::exp(-std::abs(input-mean_)/scale_)/scale_;
  }

  Real evaluateCDF(const Real input) const {
    return ((input < mean_) ? 0.5*std::exp((input-mean_)/scale_) : 
            1.-0.5*std::exp(-(input-mean_)/scale_));
  }

  Real integrateCDF(const Real input) const {
    return ((input < mean_) ? 0.5*scale_*std::exp((input-mean_)/scale_) : 
            (input-mean_)+0.5*scale_*std::exp(-(input-mean_)/scale_));
  }

  Real invertCDF(const Real input) const {
    Real sgn = ((input < 0.5) ? -1. : ((input > 0.5) ? 1. : 0.0));
    return mean_ - scale_*sgn*std::log(1.-2.*std::abs(input-0.5));
  }

  Real moment(const size_t m) const {
    if ( m == 1 ) {
      return mean_;
    }
    if ( m == 2 ) {
      return std::pow(mean_,2) + 2.*std::pow(scale_,2);
    }
    Real coeff = 0., val = 0.;
    for (size_t k = 0; k < m+1; k++) {
      if ( k%2 == 0 ) {
        coeff = compute_coeff(m,k);
        val  += coeff*std::pow(scale_,k)*std::pow(mean_,m-k);
      }
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
