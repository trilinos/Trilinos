// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ITERATIONSCALING_H
#define ROL_ITERATIONSCALING_H

/** \class ROL::IterationScaling
    \brief Provides an implementation of iteration scaled line search.
*/

#include "ROL_LineSearch.hpp"

namespace ROL { 

template<class Real>
class IterationScaling : public LineSearch<Real> {
private:
  int algo_iter_;
  ROL::Ptr<Vector<Real> > xnew_; 

public:

  virtual ~IterationScaling() {}

  // Constructor
  IterationScaling( ROL::ParameterList &parlist ) : LineSearch<Real>(parlist), algo_iter_(0) {}

  void initialize(const ROL::Vector<Real> &x, const ROL::Vector<Real> &s, const ROL::Vector<Real> &g, 
                  Objective<Real> &obj, BoundConstraint<Real> &con) {
    LineSearch<Real>::initialize(x,s,g,obj,con);
    xnew_ = x.clone();
  }

  // Run Iteration scaled line search
  void run( Real &alpha, Real &fval, int &ls_neval, int &ls_ngrad,
            const Real &gs, const Vector<Real> &s, const Vector<Real> &x, 
            Objective<Real> &obj, BoundConstraint<Real> &con ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ls_neval = 0;
    ls_ngrad = 0;
    // Get line search parameter
    algo_iter_++;
    alpha = LineSearch<Real>::getInitialAlpha(ls_neval,ls_ngrad,fval,gs,x,s,obj,con)/static_cast<Real>(algo_iter_);
    // Update iterate
    LineSearch<Real>::updateIterate(*xnew_,x,s,alpha,con);
    // Compute objective function value
    obj.update(*xnew_);
    fval = obj.value(*xnew_,tol);
    ls_neval++;
  }
};

}

#endif
