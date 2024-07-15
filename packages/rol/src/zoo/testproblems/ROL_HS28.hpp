// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 28th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS28_HPP
#define ROL_HS28_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 28th test function.
 *
 *  Exact solution x* = (0.5, -0.5, 0.5)
 *  f(x*) = 0
 */

template<class Real>
class Objective_HS28 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real two(2);
    return std::pow(x[0]+x[1],two) + std::pow(x[1]+x[2],two);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real two(2);
    g[0] = two*(x[0]+x[1]);
    g[1] = two*(x[0]+x[1]) + two*(x[1]+x[2]);
    g[2] = two*(x[1]+x[2]);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real two(2);
    hv[0] = two*v[0] + two*v[1];
    hv[1] = two*v[0] + two*two*v[1] + two*v[2];
    hv[2] = two*v[1] + two*v[2];
  } 
};

template<class Real>
class Constraint_HS28 : public StdConstraint<Real> {
public:
  Constraint_HS28(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c3(3);
    c[0] = c1*x[0]+c2*x[1]+c3*x[2]-c1;
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c1(1), c2(2), c3(3);
    jv[0] = c1*v[0]+c2*v[1]+c3*v[2];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c3(3);
    ajv[0] = c1*v[0];
    ajv[1] = c2*v[0];
    ajv[2] = c3*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }


};

template<class Real>
class getHS28 : public TestProblem<Real> {
public:
  getHS28(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS28<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 3;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(-4);
    (*xp)[1] = static_cast<Real>(1);
    (*xp)[2] = static_cast<Real>(1);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 3;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(0.5);
    (*xp)[1] = static_cast<Real>(-0.5);
    (*xp)[2] = static_cast<Real>(0.5);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return ROL::makePtr<Constraint_HS28<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    int n = 1;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
