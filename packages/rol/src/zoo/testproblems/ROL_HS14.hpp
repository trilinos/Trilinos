// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 14th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS14_HPP
#define ROL_HS14_HPP

#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 14th test function.
 *
 *  Exact solution x* = (0.5*(sqrt(7)-1), 0.25*(sqrt(7)+1))
 *  f(x*) = 9-2.875*sqrt(7)
 */

template<class Real>
class Objective_HS14 : public StdObjective<Real> {
public:

  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2);
    return std::pow(x[0]-c2,c2) + std::pow(x[1]-c1,c2);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2);
    g[0] = c2*(x[0]-c2);
    g[1] = c2*(x[1]-c1);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    hv[0] = c2*v[0];
    hv[1] = c2*v[1];
  } 
};

template<class Real>
class Constraint_HS14a : public StdConstraint<Real> {
public:
  Constraint_HS14a(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2);
    c[0] = x[0] - c2*x[1] + c1;
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c2(2);
    jv[0] = v[0] - c2*v[1];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    ajv[0] = v[0];
    ajv[1] = -c2*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }
};

template<class Real>
class Constraint_HS14b : public StdConstraint<Real> {
public:
  Constraint_HS14b(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c0(0.25), c1(1), c2(2);
    c[0] = -c0*std::pow(x[0],c2) - std::pow(x[1],c2) + c1;
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c0(0.25), c2(2);
    jv[0] = -c0*c2*x[0]*v[0] - c2*x[1]*v[1];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c0(0.25), c2(2);
    ajv[0] = -c0*c2*x[0]*v[0];
    ajv[1] = -c2*x[1]*v[0]; 
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    const Real c0(0.25), c2(2);
    ahuv[0] = -c0*c2*v[0]*u[0];
    ahuv[1] = -c2*v[1]*u[0];
  }
};

template<class Real>
class getHS14 : public TestProblem<Real> {
public:
  getHS14(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS14<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 2;
    return ROL::makePtr<StdVector<Real>>(n,2.0);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 2;
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(0.5 *(std::sqrt(7)-1));
    (*xp)[1] = static_cast<Real>(0.25*(std::sqrt(7)+1));
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return makePtr<Constraint_HS14a<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    return makePtr<StdVector<Real>>(1,0.0);
  }

  Ptr<Constraint<Real>> getInequalityConstraint(void) const {
    return makePtr<Constraint_HS14b<Real>>();
  }

  Ptr<Vector<Real>> getInequalityMultiplier(void) const {
    return makePtr<StdVector<Real>>(1,0.0);
  }

  Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
    Ptr<Vector<Real>> lb = makePtr<StdVector<Real>>(1,0.0);
    return makePtr<Bounds<Real>>(*lb,true);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
