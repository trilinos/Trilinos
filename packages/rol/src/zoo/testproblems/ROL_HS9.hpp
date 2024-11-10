// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 9th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS9_HPP
#define ROL_HS9_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 9th test function.
 *
 *  Exact solution x* = (12k-3, 16k-4), k integer
 *  f(x*) = -0.5
 */

template<class Real>
class Objective_HS9 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real pi(M_PI), c12(12), c16(16);
    return std::sin(x[0]*pi/c12)*std::cos(x[1]*pi/c16);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real pi(M_PI), c12(12), c16(16);
    g[0] =  pi/c12*std::cos(x[0]*pi/c12)*std::cos(x[1]*pi/c16);
    g[1] = -pi/c16*std::sin(x[0]*pi/c12)*std::sin(x[1]*pi/c16);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real pi(M_PI), c12(12), c16(16);
    Real h11 = -pi/c12*pi/c12*std::sin(x[0]*pi/c12)*std::cos(x[1]*pi/c16);
    Real h12 = -pi/c12*pi/c16*std::cos(x[0]*pi/c12)*std::sin(x[1]*pi/c16);
    Real h22 = -pi/c16*pi/c16*std::sin(x[0]*pi/c12)*std::cos(x[1]*pi/c16);
    hv[0] = h11*v[0] + h12*v[1];
    hv[1] = h12*v[0] + h22*v[1];
  } 
};

template<class Real>
class Constraint_HS9 : public StdConstraint<Real> {
public:
  Constraint_HS9(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c3(3), c4(4);
    c[0] = c4*x[0] - c3*x[1];
  }  

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c3(3), c4(4);
    jv[0] = c4*v[0] - c3*v[1];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c3(3), c4(4);
    ajv[0] =  c4*v[0];
    ajv[1] = -c3*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }


};

template<class Real>
class getHS9 : public TestProblem<Real> {
public:
  getHS9(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_HS9<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 2;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 2;
    // Get Solution
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(-3);
    (*xp)[1] = static_cast<Real>(-4);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    return ROL::makePtr<Constraint_HS9<Real>>();
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    // Problem dimension
    int n = 1;
    return ROL::makePtr<StdVector<Real>>(n,0.0);
  }
};



} // End ZOO Namespace
} // End ROL Namespace

#endif
