// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CUBICINTERP_H
#define ROL_CUBICINTERP_H

/** \class ROL::CubicInterp
    \brief Implements cubic interpolation back tracking line search.
*/

#include "ROL_LineSearch.hpp"

namespace ROL { 

template<class Real>
class CubicInterp : public LineSearch<Real> {
private:
  Real rho_;
  ROL::Ptr<Vector<Real> > xnew_; 

public:

  virtual ~CubicInterp() {}

  // Constructor
  CubicInterp( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist) { 
    Real half(0.5);
    rho_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").get("Backtracking Rate",half);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                   Objective<Real> &obj, BoundConstraint<Real> &con ) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
  }

  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ls_neval = 0;
    ls_ngrad = 0;
    // Get initial line search parameter
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con);
    // Update iterate
    LineSearch<Real>::updateIterate(*xnew_,x,s,alpha,con);
    // Get objective value at xnew
    Real fold = fval;
    obj.update(*xnew_);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
    // Initialize
    Real fvalp(0), alpha1(0), alpha2(0), a(0), b(0), x1(0), x2(0);
    const Real one(1), two(2), three(3), half(0.5), p1(0.1);
    bool first_iter = true;
    // Perform cubic interpolation back tracking
    while ( !LineSearch<Real>::status(LINESEARCH_CUBICINTERP,ls_neval,ls_ngrad,alpha,fold,gs,fval,x,s,obj,con) ) {
      if ( first_iter ) { // Minimize quadratic interpolate
        alpha1 = -gs*alpha*alpha/(two*(fval-fold-gs*alpha));
        first_iter = false;
      }
      else {              // Minimize cubic interpolate
        x1 = fval-fold-alpha*gs;
        x2 = fvalp-fval-alpha2*gs;
        a = (one/(alpha - alpha2))*( x1/(alpha*alpha) - x2/(alpha2*alpha2));
        b = (one/(alpha - alpha2))*(-x1*alpha2/(alpha*alpha) + x2*alpha/(alpha2*alpha2));
        if ( std::abs(a) < ROL_EPSILON<Real>() ) {
          alpha1 = -gs/(two*b);
        }
        else {
          alpha1 = (-b+sqrt(b*b-three*a*gs))/(three*a);
        }
        if ( alpha1 > half*alpha ) {
          alpha1 = half*alpha;
        }
      }
      alpha2 = alpha;
      fvalp  = fval;
      // Back track if necessary
      if ( alpha1 <= p1*alpha ) {
        alpha *= p1;
      }
      else if ( alpha1 >= half*alpha ) {
        alpha *= half;
      }
      else {
        alpha = alpha1;
      }
      // Update iterate
      LineSearch<Real>::updateIterate(*xnew_,x,s,alpha,con);
      // Get objective value at xnew
      obj.update(*xnew_);
      fval = obj.value(*xnew_,tol);
      ls_neval++;
    }
  }
};

}

#endif
