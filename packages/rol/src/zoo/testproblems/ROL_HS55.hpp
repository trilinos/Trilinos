// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 55th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS55_HPP
#define ROL_HS55_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 55th test function.
 *
 *  Exact solution x* = (0, 4/3, 5/3, 1, 2/3, 1/3)
 *  f(x*) = 19/3
 */

template<class Real>
class Objective_HS55 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c4(4);
    return x[0] + c2*x[1] + c4*x[4] + std::exp(x[0]*x[3]);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real c0(0), c1(1), c2(2), c4(4);
    Real exp03 = std::exp(x[0]*x[3]);
    g[0] = c1 + x[3]*exp03;
    g[1] = c2;
    g[2] = c0;
    g[3] = x[0]*exp03;
    g[4] = c4;
    g[5] = c0;
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real c0(0), c1(1);
    Real exp03 = std::exp(x[0]*x[3]);
    hv[0] = x[3]*x[3]*exp03*v[0] + (c1 + x[3]*x[0])*exp03*v[3];
    hv[1] = c0;
    hv[2] = c0;
    hv[3] = (c1 + x[0]*x[3])*exp03*v[0] + x[0]*x[0]*exp03*v[3];
    hv[4] = c0;
    hv[5] = c0;
  } 
};

template<class Real>
class Constraint_HS55 : public StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c3(3), c5(5), c6(6);
    c[0] = x[0] + c2*x[1] + c5*x[4] - c6;
    c[1] = x[0] + x[1] + x[2] - c3;
    c[2] = x[3] + x[4] + x[5] - c2;
    c[3] = x[0] + x[3] - c1;
    c[4] = x[1] + x[4] - c2;
    c[5] = x[2] + x[5] - c2;
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c2(2), c5(5);
    jv[0] = v[0] + c2*v[1] + c5*v[4];
    jv[1] = v[0] + v[1] + v[2];
    jv[2] = v[3] + v[4] + v[5];
    jv[3] = v[0] + v[3];
    jv[4] = v[1] + v[4];
    jv[5] = v[2] + v[5];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c2(2), c5(5);
    ajv[0] = v[0] + v[1] + v[3];
    ajv[1] = c2*v[0] + v[1] + v[4];
    ajv[2] = v[1] + v[5];
    ajv[3] = v[2] + v[3];
    ajv[4] = c5*v[0] + v[2] + v[4];
    ajv[5] = v[2] + v[5];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }
};

template<class Real>
class getHS55 : public TestProblem<Real> {
public:
  getHS55(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS55<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 6;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n);
    (*xp)[0] = static_cast<Real>(1.0);
    (*xp)[1] = static_cast<Real>(2.0);
    (*xp)[2] = static_cast<Real>(0.0);
    (*xp)[3] = static_cast<Real>(0.0);
    (*xp)[4] = static_cast<Real>(0.0);
    (*xp)[5] = static_cast<Real>(2.0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 6;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n);
    (*xp)[0] = static_cast<Real>(0.0);
    (*xp)[1] = static_cast<Real>(4.0/3.0);
    (*xp)[2] = static_cast<Real>(5.0/3.0);
    (*xp)[3] = static_cast<Real>(1.0);
    (*xp)[4] = static_cast<Real>(2.0/3.0);
    (*xp)[5] = static_cast<Real>(1.0/3.0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return ROL::makePtr<Constraint_HS55<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    int n = 6;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    int n = 6;
    Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(n,0.0);
    ROL::Ptr<std::vector<Real>> up = ROL::makePtr<std::vector<Real>>(n,ROL_INF<Real>());
    (*up)[0] = static_cast<Real>(1);
    (*up)[3] = static_cast<Real>(1);
    Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(up);
    return makePtr<Bounds<Real>>(l,u);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
