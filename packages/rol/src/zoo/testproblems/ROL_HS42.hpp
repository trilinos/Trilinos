// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 42th test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS42_HPP
#define ROL_HS42_HPP

#include "ROL_StdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 42th test function.
 *
 *  Exact solution x* = (2, 2, 0.6*sqrt(2), 0.8*sqrt(2))
 *  f(x*) = 28-10*sqrt(2)
 */

template<class Real>
class Objective_HS42 : public StdObjective<Real> {
public:

  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c3(3), c4(4);
    return std::pow(x[0]-c1,c2) + std::pow(x[1]-c2,c2)
         + std::pow(x[2]-c3,c2) + std::pow(x[3]-c4,c2);
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real c1(1), c2(2), c3(3), c4(4);
    g[0] = c2*(x[0]-c1);
    g[1] = c2*(x[1]-c2);
    g[2] = c2*(x[2]-c3);
    g[3] = c2*(x[3]-c4);
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    hv[0] = c2*v[0];
    hv[1] = c2*v[1];
    hv[2] = c2*v[2];
    hv[3] = c2*v[3];
  }
};

template<class Real>
class Constraint_HS42a : public StdConstraint<Real> {
public:
  Constraint_HS42a(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    c[0] = x[0] - c2;
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    jv[0] = v[0];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    ajv[0] = v[0];
    ajv[1] = static_cast<Real>(0);
    ajv[2] = static_cast<Real>(0);
    ajv[3] = static_cast<Real>(0);
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }
};

template<class Real>
class Constraint_HS42b : public StdConstraint<Real> {
public:
  Constraint_HS42b(void) {}
 
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    c[0] = std::pow(x[2],c2) + std::pow(x[3],c2) - c2;
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c2(2);
    jv[0] = c2*x[2]*v[2] + c2*x[3]*v[3];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    ajv[0] = static_cast<Real>(0);
    ajv[1] = static_cast<Real>(0);
    ajv[2] = c2*x[2]*v[0];
    ajv[3] = c2*x[3]*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    const Real c2(2);
    ahuv[0] = static_cast<Real>(0);
    ahuv[1] = static_cast<Real>(0);
    ahuv[2] = c2*v[2]*u[0];
    ahuv[3] = c2*v[3]*u[0];
  }
};

template<class Real>
class getHS42 : public TestProblem<Real> {
public:
  getHS42(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS42<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 4;
    return ROL::makePtr<StdVector<Real>>(n,1.0);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 4;
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(2);
    (*xp)[1] = static_cast<Real>(2);
    (*xp)[2] = static_cast<Real>(0.6*std::sqrt(2));
    (*xp)[3] = static_cast<Real>(0.8*std::sqrt(2));
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    std::vector<Ptr<Constraint<Real>>> cvec(2);
    cvec[0] = makePtr<Constraint_HS42a<Real>>();
    cvec[1] = makePtr<Constraint_HS42b<Real>>();
    return ROL::makePtr<Constraint_Partitioned<Real>>(cvec);
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    std::vector<Ptr<Vector<Real>>> lvec(2);
    lvec[0] = makePtr<StdVector<Real>>(makePtr<std::vector<Real>>(1,0.0));
    lvec[1] = makePtr<StdVector<Real>>(makePtr<std::vector<Real>>(1,0.0));
    return ROL::makePtr<PartitionedVector<Real>>(lvec);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
