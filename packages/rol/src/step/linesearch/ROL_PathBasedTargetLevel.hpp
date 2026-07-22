// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PATHBASEDTARGETLEVEL_H
#define ROL_PATHBASEDTARGETLEVEL_H

/** \class ROL::PathBasedTargetLevel
    \brief Provides an implementation of path-based target leve line search.
*/

#include "ROL_LineSearch.hpp"

namespace ROL { 

template<class Real>
class PathBasedTargetLevel : public LineSearch<Real> {
private:
  ROL::Ptr<Vector<Real> > xnew_; 

  Real min_value_;
  Real rec_value_;
  Real target_;
  Real delta_;
  Real sigma_;
  Real bound_;

public:

  virtual ~PathBasedTargetLevel() {}

  // Constructor
  PathBasedTargetLevel( ROL::ParameterList &parlist ) 
    : LineSearch<Real>(parlist), min_value_(ROL::ROL_OVERFLOW<Real>()),
      rec_value_(ROL::ROL_OVERFLOW<Real>()),  target_(0.0), sigma_(0.0) {
    Real p1(0.1), one(1);
    delta_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").sublist("Path-Based Target Level").get("Target Relaxation Parameter",p1);
    bound_ = parlist.sublist("Step").sublist("Line Search").sublist("Line-Search Method").sublist("Path-Based Target Level").get("Upper Bound on Path Length",one);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g,
                  Objective<Real> &obj, BoundConstraint<Real> &con) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
  }

  // Run Iteration scaled line search
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>()), zero(0), half(0.5);
    ls_neval = 0;
    ls_ngrad = 0;
    // Update target objective value
    if ( fval < min_value_ ) {
      min_value_ = fval;
    }
    target_ = rec_value_ - half*delta_;
    if ( fval < target_ ) {
      rec_value_ = min_value_; 
      sigma_ = zero;
    }
    else {
      if ( sigma_ > bound_ ) {
        rec_value_ = min_value_;
        sigma_ = zero;
        delta_ *= half;
      }
    }
    target_ = rec_value_ - delta_;
    // Get line-search parameter
    alpha = (fval - target_)/std::abs(gs);
    // Update iterate
    LineSearch<Real>::updateIterate(*xnew_,x,s,alpha,con);
    // Compute objective function value
    obj.update(*xnew_);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
    // Update sigma 
    sigma_ += alpha*std::sqrt(std::abs(gs));
  }
};

}

#endif
