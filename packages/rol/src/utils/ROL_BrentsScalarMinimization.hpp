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

#ifndef ROL_BRENTSSCALARMINIMIZATION_H
#define ROL_BRENTSSCALARMINIMIZATION_H

/** \class ROL::BrentsScalarMinimization
    \brief Implements Brent's method to minimize a scalar function on a
           bounded interval.
*/

#include "ROL_ScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class BrentsScalarMinimization : public ScalarMinimization<Real> {
private:
  Real tol_;
  int niter_;

public:
  // Constructor
  BrentsScalarMinimization( ROL::ParameterList &parlist ) {
    ROL::ParameterList &list = parlist.sublist("Scalar Minimization").sublist("Brent's");
    tol_   = list.get("Tolerance",1.e-10);
    niter_ = list.get("Iteration Limit",1000);
  }

  void run(Real &fx, Real &x, int &nfval, int &ngrad,
           ScalarFunction<Real> &f, const Real A, const Real B,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real zero(0), half(0.5), one(1), two(2), three(3), five(5);
    nfval = 0; ngrad = 0;
    // ---> Set algorithmic constants
    const Real c   = half*(three - std::sqrt(five));
    const Real eps = std::sqrt(ROL_EPSILON<Real>());
    // ---> Set end points and initial guess
    Real a = A, b = B;
    x = a + c*(b-a);
    // ---> Evaluate function
    fx = f.value(x);
    nfval++;
    // ---> Initialize algorithm storage
    Real v = x, w = v, u = zero, fu = zero;
    Real p = zero, q = zero, r = zero, d = zero, e = zero;
    Real fv = fx, fw = fx, tol = zero, t2 = zero, m = zero, gx = ROL_INF<Real>();
    bool deriv = false;
    for (int i = 0; i < niter_; i++) {
      m = half*(a+b);
      tol = eps*std::abs(x) + tol_; t2 = two*tol;
      // Check stopping criterion
      if (std::abs(x-m) <= t2 - half*(b-a) || test.check(x,fx,gx,nfval,ngrad,deriv)) {
        break;
      }
      p = zero; q = zero; r = zero;
      if ( std::abs(e) > tol ) {
        // Fit parabola
        r = (x-w)*(fx-fv);     q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r; q = two*(q-r);
        if ( q > zero ) {
          p *= -one;
        }
        q = std::abs(q);
        r = e; e = d;
      }
      if ( std::abs(p) < std::abs(half*q*r) && p > q*(a-x) && p < q*(b-x) ) {
        // A parabolic interpolation step
        d = p/q; u = x + d;
        // f must not be evaluated too close to a or b
        if ( (u - a) < t2 || (b - u) < t2 ) {
          d = (x < m) ? tol : -tol;
        }
      }
      else {
        // A golden section step
        e = ((x < m) ? b : a) - x; d = c*e;
      }
      // f must not be evaluated too close to x
      u  = x + ((std::abs(d) >= tol) ? d : ((d > zero) ? tol : -tol));
      fu = f.value(u);
      nfval++;
      // Update a, b, v, w, and x
      if ( fu <= fx ) {
        if ( u < x ) {
          b = x;
        }
        else {
          a = x;
        }
        v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
      }
      else {
        if ( u < x ) {
          a = u;
        }
        else {
          b = u;
        }
        if ( fu <= fw || w == x ) {
          v = w; fv = fw; w = u; fw = fu;
        }
        else if ( fu <= fv || v == x || v == w ) {
          v = u; fv = fu;
        }
      }
    }
  }
};

}

#endif
