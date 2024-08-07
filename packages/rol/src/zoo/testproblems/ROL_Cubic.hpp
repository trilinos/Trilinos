// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for a cubic test problem.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_CUBIC_HPP
#define ROL_CUBIC_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

  template<class Real>
  class Objective_Cubic : public StdObjective<Real> {
  public:
    Objective_Cubic() {}

    Real value( const std::vector<Real> &x, Real &tol ) {
      return std::pow(x[0],3)+std::pow(x[1],3);
    }

    void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
      const Real three(3);
      g[0] = three*std::pow(x[0],2);
      g[1] = three*std::pow(x[1],2);
    }
#if USE_HESSVEC
    void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real six(6);
      hv[0] = six*x[0]*v[0];
      hv[1] = six*x[1]*v[1];
    }
#endif
    void invHessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real six(6);
      hv[0] = v[0]/(six*x[0]);
      hv[1] = v[1]/(six*x[1]);
    }
  };

  template<class Real>
  class Constraint_Cubic : public StdConstraint<Real> {
  public:
    Constraint_Cubic() {}

    void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
      c[0] = std::pow(x[0],3) + x[1];
    }

    void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real three(3);
      jv[0] = three*std::pow(x[0],2)*v[0] + v[1];
    }

    void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real three(3);
      ajv[0] = three*std::pow(x[0],2)*v[0];
      ajv[1] = v[0];
    }
#if USE_HESSVEC
    void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real zero(0), six(6);
      ahuv[0] = six*x[0]*u[0]*v[0];
      ahuv[1] = zero;
    }
#endif
  };

  template<class Real>
  class getCubic : public TestProblem<Real> {
  private:
    const int type_;

  public:
    getCubic(int type = 0) : type_(type) {}

    Ptr<Objective<Real>> getObjective(void) const {
      return makePtr<Objective_Cubic<Real>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>( 1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(-0.9)); 
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>( 1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(-1.0));
      if (type_ == 1) {
        const Real one(1), /*two(2),*/ three(3), six(6);
        Real x = -one/std::pow(three,one/six);
        Real y = -std::pow(x,3);
        (*xp)[0] = x;
        (*xp)[1] = y;
      }
      if (type_ == 2) {
        // This solution is only approximate
        (*xp)[0] = static_cast<Real>(-0.8374930678347255);
        (*xp)[1] = static_cast<Real>( 0.5774131462277658);
      }
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
      int n = 2;
      Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(n,-1.0);
      Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(n, 1.0);
      return makePtr<Bounds<Real>>(l,u);
    }

    Ptr<Constraint<Real>> getEqualityConstraint(void) const {
      if (type_ == 1) {
        return makePtr<Constraint_Cubic<Real>>();
      }
      return nullPtr;
    }

    Ptr<Vector<Real>> getEqualityMultiplier(void) const {
      if (type_ == 1) {
        return makePtr<StdVector<Real>>(1,0.0);
      }
      return nullPtr;
    }

    Ptr<Constraint<Real>> getInequalityConstraint(void) const {
      if (type_ == 2) {
        return makePtr<Constraint_Cubic<Real>>();
      }
      return nullPtr;
    }

    Ptr<Vector<Real>> getInequalityMultiplier(void) const {
      if (type_ == 2) {
        return makePtr<StdVector<Real>>(1,0.0);
      }
      return nullPtr;
    }

    Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
      if (type_ == 2) {
        Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(1,-0.01);
        Ptr<Vector<Real>> u = makePtr<StdVector<Real>>(1, 0.01);
        return makePtr<Bounds<Real>>(l,u);
      }
      return nullPtr;
    }
  };

}// End ZOO Namespace
}// End ROL Namespace

#endif
