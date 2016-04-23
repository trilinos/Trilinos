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

#ifndef ROL_MEANVARIANCEQUAD_HPP
#define ROL_MEANVARIANCEQUAD_HPP

#include "ROL_ExpectationQuad.hpp"

namespace ROL {

template<class Real>
class MeanVarianceQuadrangle : public ExpectationQuad<Real> {
private:
  Real coeff_;

  void checkInputs(void) const {
    Real zero(0);
    TEUCHOS_TEST_FOR_EXCEPTION((coeff_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::MeanVarianceQuadrangle): Coefficient must be positive!");
  }

public:

  MeanVarianceQuadrangle(const Real coeff = 1)
    : ExpectationQuad<Real>(), coeff_(coeff) {
    checkInputs();
  }

  MeanVarianceQuadrangle(Teuchos::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mean-Variance Quadrangle");
    coeff_ = list.get<Real>("Coefficient");
    checkInputs();
  }

  Real error(Real x, int deriv = 0) {
    Real err(0), two(2);
    if (deriv==0) {
      err = coeff_*x*x;
    }
    else if (deriv==1) {
      err = two*coeff_*x;
    }
    else {
      err = two*coeff_;
    }
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    Real reg = error(x,deriv) + X;
    return reg;
  }

};

}
#endif
