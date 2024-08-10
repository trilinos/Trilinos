// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 52nd test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS52_HPP
#define ROL_HS52_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 52nd test function.
 *
 *  Exact solution x* = (1, 1, 1, 1, 1)
 *  f(x*) = 0
 */

template<class Real>
class Objective_HS52 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c4(4);
    return std::pow(c4*x[0]-x[1],c2) + std::pow(x[1]+x[2]-c2,c2)
         + std::pow(x[3]-c1,c2) + std::pow(x[4]-c1,c2);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c4(4);
    g[0] = c2*c4*(c4*x[0]-x[1]);
    g[1] = c2*(x[1]-c4*x[0]) + c2*(x[1]+x[2]-c2);
    g[2] = c2*(x[1]+x[2]-c2);
    g[3] = c2*(x[3]-c1);
    g[4] = c2*(x[4]-c1);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c4(4);
    hv[0] = c2*c4*c4*v[0] - c2*c4*v[1];
    hv[1] = (c2+c2)*v[1] - c2*c4*v[0] + c2*v[2];
    hv[2] = c2*v[2] + c2*v[1];
    hv[3] = c2*v[3];
    hv[4] = c2*v[4];
  } 
};

template<class Real>
class Constraint_HS52 : public StdConstraint<Real> {
public:
  Constraint_HS52(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c3(3);
    c[0] = x[0]+c3*x[1];
    c[1] = x[2]+x[3]-c2*x[4];
    c[2] = x[1]-x[4];
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c2(2), c3(3);
    jv[0] = v[0]+c3*v[1];
    jv[1] = v[2]+v[3]-c2*v[4];
    jv[2] = v[1]-v[4];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c3(3);
    ajv[0] = v[0];
    ajv[1] = c3*v[0] + v[2];
    ajv[2] = v[1];
    ajv[3] = v[1];
    ajv[4] = -c2*v[1] - v[2];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }


};

template<class Real>
class getHS52 : public TestProblem<Real> {
public:
  getHS52(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS52<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 5;
    return ROL::makePtr<StdVector<Real>>(n,2.0);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 5;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(-33.0/349.0);
    (*xp)[1] = static_cast<Real>(11.0/349.0);
    (*xp)[2] = static_cast<Real>(180.0/349.0);
    (*xp)[3] = static_cast<Real>(-158.0/349.0);
    (*xp)[4] = static_cast<Real>(11.0/349.0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return ROL::makePtr<Constraint_HS52<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    int n = 3;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
