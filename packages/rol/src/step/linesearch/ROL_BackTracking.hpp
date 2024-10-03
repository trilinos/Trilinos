// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BACKTRACKING_H
#define ROL_BACKTRACKING_H

#include "ROL_LineSearch.hpp"

/** \class ROL::BackTracking
    \brief Implements a simple back tracking line search.
*/

namespace ROL { 

template<class Real>
class BackTracking : public LineSearch<Real> {
private:
  Real rho_;
  ROL::Ptr<Vector<Real> > xnew_; 

public:

  virtual ~BackTracking() {}

  // Constructor
  BackTracking( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist) {
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
    // Perform backtracking
    while ( !LineSearch<Real>::status(LINESEARCH_BACKTRACKING,ls_neval,ls_ngrad,alpha,fold,gs,fval,*xnew_,s,obj,con) ) {
      alpha *= rho_;
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
