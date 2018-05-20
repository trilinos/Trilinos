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

#ifndef ROL_GOLDENSECTIONSCALARMINIMIZATION_H
#define ROL_GOLDENSECTIONSCALARMINIMIZATION_H

/** \class ROL::GoldenSectionScalarMinimization
    \brief Implements the golden section method to minimize a scalar function on a
           bounded interval.
*/

#include "ROL_ScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class GoldenSectionScalarMinimization : public ScalarMinimization<Real> {
private:
  Real tol_;
  int niter_;

public:
  // Constructor
  GoldenSectionScalarMinimization( ROL::ParameterList &parlist ) {
    ROL::ParameterList &list = parlist.sublist("Scalar Minimization").sublist("Golden Section");
    tol_   = list.get("Tolerance",1.e-10);
    niter_ = list.get("Iteration Limit",1000);
  }

  void run(Real &fx, Real &x, int &nfval, int &ngrad,
           ScalarFunction<Real> &f,
           const Real A, const Real B,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real one(1), two(2), five(5);
    nfval = 0; ngrad = 0;
    // Reciprocal of golden ratio
    const Real c = two/(one + std::sqrt(five));
    // Compute value f(a), f(b), f(u), and f(v)
    Real a = A,               fa = f.value(a); nfval++;
    Real b = B,               fb = f.value(b); nfval++;
    Real u = c*a + (one-c)*b, fu = f.value(u); nfval++;
    Real v = (one-c)*a + c*b, fv = f.value(v); nfval++;
    Real gx = ROL_INF<Real>();
    bool deriv = false;
    // Find minimum of all function evaluations
    if ( fa <= fu && fa <= fv && fa <= fb ) {
      x = a; fx = fa;
    }
    else if ( fu <= fa && fu <= fv && fu <= fb ) {
      x = u; fx = fu;
    }
    else if ( fv <= fa && fv <= fu && fv <= fb ) {
      x = v; fx = fv;
    }
    else {
      x = b; fx = fb;
    }
    // Run Golden Section
    for ( int i = 0; i < niter_; i++ ) {
      if ( std::abs(b - a) < tol_ || test.check(x,fx,gx,nfval,ngrad,deriv) ) {
        break;
      }
      if ( fu > fv ) {
        a = u;               fa = fu;
        u = v;               fu = fv;
        v = (one-c)*a + c*b; fv = f.value(v); nfval++;
      }
      else {
        b = v;               fb = fv;
        v = u;               fv = fu;
        u = c*a + (one-c)*b; fu = f.value(u); nfval++;
      }
      if ( fa <= fu && fa <= fv && fa <= fb ) {
        x = a; fx = fa;
      }
      else if ( fu <= fa && fu <= fv && fu <= fb ) {
        x = u; fx = fu;
      }
      else if ( fv <= fa && fv <= fu && fv <= fb ) {
        x = v; fx = fv;
      }
      else {
        x = b; fx = fb;
      }
    }
  }
};

}

#endif
