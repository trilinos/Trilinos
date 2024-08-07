// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for the cylinder head test problem.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_CANTILEVER_HPP
#define ROL_CANTILEVER_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

  template<class Real>
  class Objective_Cantilever : public StdObjective<Real> {
  public:
    Objective_Cantilever() {}

    Real value( const std::vector<Real> &x, Real &tol ) {
      return x[0]*x[1];
    }

    void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
      g[0] = x[1];
      g[1] = x[0];
    }
#if USE_HESSVEC
    void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      hv[0] = v[1];
      hv[1] = v[0];
    }
#endif
  };

  template<class Real>
  class Constraint_Cantilever : public StdConstraint<Real> {
  private:
    Real stress(const Real w, const Real t, const int deriv = 0, const int comp1 = 0, const int comp2 = 0) const {
      const Real scale(600), X(500), Y(1000);
      Real val(0);
      if (deriv == 0) {
        val = scale*(Y/(w*t*t) + X/(w*w*t));
      }
      else if (deriv == 1) {
        if (comp1 == 0) {
          const Real two(2);
          val = scale*(-Y/(w*w*t*t) - two*X/(w*w*w*t));
        }
        else if (comp1 == 1) {
          const Real two(2);
          val = scale*(-two*Y/(w*t*t*t) - X/(w*w*t*t));
        }
      }
      else if (deriv == 2) {
        if (comp1 == 0 && comp2 == 0) {
          const Real two(2), six(6);
          val = scale*(two*Y/(w*w*w*t*t) + six*X/(w*w*w*w*t));
        }
        else if (comp1 == 1 && comp2 == 1) {
          const Real two(2), six(6);
          val = scale*(six*Y/(w*t*t*t*t) + two*X/(w*w*t*t*t));
        }
        else if (comp1 == 0 && comp2 == 1) {
          const Real two(2);
          val = scale*two*(Y/(w*w*t*t*t) + X/(w*w*w*t*t));
        }
        else if (comp1 == 1 && comp2 == 0) {
          const Real two(2);
          val = scale*two*(Y/(w*w*t*t*t) + X/(w*w*w*t*t));
        }
      }
      return val;
    }

    Real displacement(const Real w, const Real t, const int deriv = 0, const int comp1 = 0, const int comp2 = 0) const {
      const Real four(4), L(100), E(2.9e7), X(500), Y(1000);
      const Real C = four*std::pow(L,3)/E;
      Real arg1 = std::pow(Y/(t*t),2), arg2 = std::pow(X/(w*w),2);
      Real mag  = std::sqrt(arg1 + arg2);
      Real val(0);
      if (deriv == 0) {
        val = C/(w*t)*mag;
      }
      else if (deriv == 1) {
        if (comp1 == 0) {
          const Real three(3);
          val = -C * (three * std::pow(X*t*t,2) + std::pow(Y*w*w,2))
                    / (std::pow(w,6)*std::pow(t,5)*mag);
        }
        else if (comp1 == 1) {
          const Real three(3);
          val = -C * (std::pow(X*t*t,2) + three*std::pow(Y*w*w,2))
                    / (std::pow(w,5)*std::pow(t,6)*mag);
        }
      }
      else if (deriv == 2) {
        if (comp1 == 0 && comp2 == 0) {
          const Real two(2), six(6), nine(9);
          val = C * two * mag * (std::pow(Y*w*w,4) + nine*std::pow(Y*X*w*w*t*t,2) + six*std::pow(X*t*t,4))
                   / (std::pow(w,3)*t*std::pow(std::pow(Y*w*w,2)+std::pow(X*t*t,2),2));
        }
        else if (comp1 == 1 && comp2 == 1) {
          const Real two(2), six(6), nine(9);
          val = C * two * mag * (six*std::pow(Y*w*w,4) + nine*std::pow(Y*X*w*w*t*t,2) + std::pow(X*t*t,4))
                   / (std::pow(t,3)*w*std::pow(std::pow(Y*w*w,2)+std::pow(X*t*t,2),2));
        }
        else if (comp1 == 0 && comp2 == 1) {
          const Real two(2), three(3);
          val = C * (three*std::pow(X*t*t,4) + two*std::pow(X*Y*t*t*w*w,2) + three*std::pow(Y*w*w,4))
                    / (std::pow(t*w,6)*mag*(std::pow(X*t*t,2) + std::pow(Y*w*w,2)));
        }
        else if (comp1 == 1 && comp2 == 0) {
          const Real two(2), three(3);
          val = C * (three*std::pow(X*t*t,4) + two*std::pow(X*Y*t*t*w*w,2) + three*std::pow(Y*w*w,4))
                    / (std::pow(t*w,6)*mag*(std::pow(X*t*t,2) + std::pow(Y*w*w,2)));
        }
      }
      return val;
    }
  public:
    Constraint_Cantilever() {}

    void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
      const Real R(40000), D(2.2535), one(1);
      Real s = stress(x[0],x[1],0)/R;
      Real d = displacement(x[0],x[1],0)/D;
      c[0] = s - one;
      c[1] = d - one;
    }

    void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real R(40000), D(2.2535);
      Real s0 = stress(x[0],x[1],1,0)/R, s1 = stress(x[0],x[1],1,1)/R;
      Real d0 = displacement(x[0],x[1],1,0)/D, d1 = displacement(x[0],x[1],1,1)/D;
      jv[0] = s0*v[0] + s1*v[1];
      jv[1] = d0*v[0] + d1*v[1];
    }

    void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real R(40000), D(2.2535);
      Real s0 = stress(x[0],x[1],1,0)/R, s1 = stress(x[0],x[1],1,1)/R;
      Real d0 = displacement(x[0],x[1],1,0)/D, d1 = displacement(x[0],x[1],1,1)/D;
      ajv[0] = s0*v[0] + d0*v[1];
      ajv[1] = s1*v[0] + d1*v[1];
    }
#if USE_HESSVEC
    void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real R(40000), D(2.2535);
      Real s00 = stress(x[0],x[1],2,0,0)/R, s01 = stress(x[0],x[1],2,0,1)/R;
      Real s10 = stress(x[0],x[1],2,1,0)/R, s11 = stress(x[0],x[1],2,1,1)/R;
      Real d00 = displacement(x[0],x[1],2,0,0)/D, d01 = displacement(x[0],x[1],2,0,1)/D;
      Real d10 = displacement(x[0],x[1],2,1,0)/D, d11 = displacement(x[0],x[1],2,1,1)/D;
      ahuv[0] = (s00*u[0] + d00*u[1])*v[0] + (s01*u[0] + d01*u[1])*v[1];
      ahuv[1] = (s10*u[0] + d10*u[1])*v[0] + (s11*u[0] + d11*u[1])*v[1];
    }
#endif
  };

  template<class Real>
  class getCantilever : public TestProblem<Real> {
  public:
    getCantilever() {}

    Ptr<Objective<Real>> getObjective(void) const {
      return makePtr<Objective_Cantilever<Real>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0)); 
      (*xp)[0] = static_cast<Real>(2.0);
      (*xp)[1] = static_cast<Real>(2.0);
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      (*xp)[0] = static_cast<Real>(2.3520341271);
      (*xp)[1] = static_cast<Real>(3.3262784077);
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> lp    = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> up    = makePtr<std::vector<Real>>(n,static_cast<Real>(4.0));
      Ptr<Vector<Real>> l = makePtr<PrimalScaledStdVector<Real>>(lp,scale);
      Ptr<Vector<Real>> u = makePtr<PrimalScaledStdVector<Real>>(up,scale);
      return makePtr<Bounds<Real>>(l,u);
    }

    Ptr<Constraint<Real>> getInequalityConstraint(void) const {
      return makePtr<Constraint_Cantilever<Real>>();
    }

    Ptr<Vector<Real>> getInequalityMultiplier(void) const {
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(2,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> lp    = makePtr<std::vector<Real>>(2,static_cast<Real>(0.0)); 
      return makePtr<DualScaledStdVector<Real>>(lp,scale);
    }

    Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(2,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> up    = makePtr<std::vector<Real>>(2,static_cast<Real>(0.0)); 
      Ptr<Vector<Real>> u = makePtr<DualScaledStdVector<Real>>(up,scale);
      return makePtr<Bounds<Real>>(*u,false);
    }
  };

}// End ZOO Namespace
}// End ROL Namespace

#endif
