// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 48th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS48_HPP
#define ROL_HS48_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 48th test function.
 *
 *  Exact solution x* = (1, 1, 1, 1, 1)
 *  f(x*) = 0
 */

template<class Real>
class Objective_HS48 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real one(1), two(2);
    return std::pow(x[0]-one,two) + std::pow(x[1]-x[2],two) + std::pow(x[3]-x[4],two);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real one(1), two(2);
    g[0] = two*(x[0]-one);
    g[1] = two*(x[1]-x[2]);
    g[2] = two*(x[2]-x[1]);
    g[3] = two*(x[3]-x[4]);
    g[4] = two*(x[4]-x[3]);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real two(2);
    hv[0] = two*v[0];
    hv[1] = two*v[1] - two*v[2];
    hv[2] = two*v[2] - two*v[1];
    hv[3] = two*v[3] - two*v[4];
    hv[4] = two*v[4] - two*v[3];
  } 
};

template<class Real>
class Constraint_HS48 : public StdConstraint<Real> {
public:
  Constraint_HS48(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c3(3), c5(5);
    c[0] = x[0]+x[1]+x[2]+x[3]+x[4]-c5;
    c[1] = x[2]-c2*(x[3]+x[4])+c3;
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c2(2);
    jv[0] = v[0]+v[1]+v[2]+v[3]+v[4];
    jv[1] = v[2]-c2*(v[3]+v[4]);
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    ajv[0] = v[0];
    ajv[1] = v[0];
    ajv[2] = v[0] + v[1];
    ajv[3] = v[0] - c2*v[1];
    ajv[4] = v[0] - c2*v[1];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }


};

template<class Real>
class getHS48 : public TestProblem<Real> {
public:
  getHS48(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS48<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 5;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(3);
    (*xp)[1] = static_cast<Real>(5);
    (*xp)[2] = static_cast<Real>(-3);
    (*xp)[3] = static_cast<Real>(2);
    (*xp)[4] = static_cast<Real>(-2);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 5;
    return ROL::makePtr<StdVector<Real>>(n,1.0);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return ROL::makePtr<Constraint_HS48<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    int n = 2;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
