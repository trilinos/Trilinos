// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for sum of squares function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_SUMOFSQUARES_HPP
#define ROL_SUMOFSQUARES_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Sum of squares function. 
 */
template<class Real>
class Objective_SumOfSquares : public Objective<Real> {
public:
  Real value( const Vector<Real> &x, Real &tol ) {
    return x.dot(x);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    g.set(x);
    g.scale(2.0);
  }

#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.set(v);
    hv.scale(2.0);
  }
#endif

  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.set(v);
    hv.scale(0.5);
  }
};

template<class Real>
class getSumOfSquares : public TestProblem<Real> {
public:
  getSumOfSquares(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_SumOfSquares<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 100;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,1.0);
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 100;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
