// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BISECTIONSCALARMINIMIZATION_H
#define ROL_BISECTIONSCALARMINIMIZATION_H

/** \class ROL::BisectionScalarMinimization
    \brief Implements the bisection method to minimize a scalar function on a
           bounded interval.
*/

#include "ROL_ScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class BisectionScalarMinimization : public ScalarMinimization<Real> {
private:
  Real tol_;
  int niter_;

public:
  // Constructor
  BisectionScalarMinimization( ROL::ParameterList &parlist ) {
    ROL::ParameterList &list = parlist.sublist("Scalar Minimization").sublist("Bisection");
    tol_   = list.get("Tolerance",1.e-10);
    niter_ = list.get("Iteration Limit",1000);
  }

  void run(Real &fx, Real &x, int &nfval, int &ngrad,
           ScalarFunction<Real> &f, const Real A, const Real B,
           ScalarMinimizationStatusTest<Real> &test) const {
    Real zero(0), half(0.5);
    nfval = 0; ngrad = 0;
    // Compute value f(A) and f(B)
    Real a = A,         fa = f.value(a); nfval++;
    Real b = B,         fb = f.value(b); nfval++;
    Real m = half*(A+B), fm = f.value(m); nfval++;
    Real u = zero, fu = zero, v = zero, fv = zero, gx = ROL_INF<Real>();
    bool deriv = false;
    // Get minimum of all evaluations
    if ( fa <= fm && fa <= fb ) {
      x = a; fx = fa;
    }
    else if ( fm <= fa && fm <= fb ) {
      x = m; fx = fm;
    }
    else {
      x = b; fx = fb;
    }
    // Run bisection
    for ( int i = 0; i < niter_; i++ ) {
      if ( std::abs(b - a) < tol_ || test.check(x,fx,gx,nfval,ngrad,deriv) ) {
        break;
      }
      u = half*(a+m); fu = f.value(u); nfval++;
      v = half*(m+b); fv = f.value(v); nfval++;

      if (    ( (fa <= fb) && (fa <= fu) && (fa <= fv) && (fa <= fm) )
           || ( (fu <= fb) && (fu <= fa) && (fu <= fv) && (fu <= fm) ) ) {
        if ( fa < fu ) {
          x = a; fx = fa;
        }
        else {
          x = u; fx = fu;
        }
        b = m; fb = fm;
        m = u; fm = fu;
      }
      else if ( ( (fm <= fb) && (fm <= fa) && (fm <= fu) && (fm <= fv) ) ) {
        x = m; fx = fm;
        a = u; fa = fu;
        b = v; fb = fv;
      }
      else if (    ( (fv <= fb) && (fv <= fa) && (fv <= fu) && (fv <= fm) )
                || ( (fb <= fa) && (fb <= fu) && (fb <= fv) && (fb <= fm) ) ) {
        if ( fb < fv ) {
          x = b; fx = fb;
        }
        else {
          x = v; fx = fv;
        }
        a = m; fa = fm;
        m = v; fm = fv;
      }
    }
  }
};

}

#endif
