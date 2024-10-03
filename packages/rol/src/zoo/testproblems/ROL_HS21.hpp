// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
 *  \brief Contains definitions for W. Hock and K. Schittkowski 21th test problem
 *         which contains bound and inequality constraints.
 */

#ifndef ROL_HS21_HPP
#define ROL_HS21_HPP

#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"

namespace ROL {
namespace ZOO {

template<class Real>
class Objective_HS21 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c1(0.1), c2(100);
    return c1*x[0]*x[0] + x[1]*x[1] - c2;
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real two(2), c1(0.1);
    g[0] = two*c1*x[0];
    g[1] = two*x[1];
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real two(2), c1(0.1);
    hv[0] = two*c1*v[0];
    hv[1] = two*v[1];
  }
}; // class Objective_HS21

template<class Real>
class Constraint_HS21 : public StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(10);
    c[0] = c1*x[0] - x[1] - c1;
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v,
                      const std::vector<Real> &x, Real &tol ) {
    const Real c1(10);
    jv[0] = c1*v[0] - v[1];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c1(10);
    ajv[0] = c1*v[0];
    ajv[1] = -v[0];
  }
 
  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                            const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }
}; // class Constraint_HS21

template<class Real>
class getHS21 : public TestProblem<Real> {
public:
  getHS21(void) {}

  Ptr<Objective<Real> > getObjective( void ) const {
    return makePtr<Objective_HS21<Real>>();
  }

  Ptr<Constraint<Real> > getInequalityConstraint( void ) const {
    return makePtr<Constraint_HS21<Real>>();
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint( void ) const {
    Ptr<std::vector<Real>> lp = makePtr<std::vector<Real>>(2,2.0);
    (*lp)[1] = static_cast<Real>(-50);

    Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(2,50.0);

    return makePtr<Bounds<Real>>(l,u);
  }

  Ptr<Vector<Real>> getInitialGuess( void ) const {
    return makePtr<StdVector<Real>>(2,-1.0);
  }

  Ptr<Vector<Real>> getSolution( const int i = 0 ) const {
    Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(2,0.0);
    (*xp)[0] = static_cast<Real>(2);
    return makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real>> getInequalityMultiplier( void ) const {
    return makePtr<StdVector<Real>>(1,0.0);
  }

  Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
    Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(1,0.0);
    Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(1,ROL_INF<Real>());
    return makePtr<Bounds<Real>>(l,u);
  }
};

} // namespace ZOO
} // namespace ROL

#endif // ROL_HS21_HPP
