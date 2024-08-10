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

#ifndef ROL_CYLINDERHEAD_HPP
#define ROL_CYLINDERHEAD_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

  template<class Real>
  class Objective_CylinderHead : public StdObjective<Real> {
  private:
    Real warranty(const Real flatness, const int deriv = 0) const {
      const Real w0(1e5), w1(1.5e4), w2(4), zero(0);
      return (deriv==0 ? w0 + w1*(w2 - flatness)
            :(deriv==1 ? -w1 : zero));
    }

    Real horsepower(const Real dintake, const int deriv = 0) const {
      const Real d0(250), d1(200), d2(1), d3(1.833), zero(0);
      return (deriv==0 ? d0 + d1*(dintake/d3 - d2)
            :(deriv==1 ? d1/d3 : zero));
    }

  public:
    Objective_CylinderHead() {}

    Real value( const std::vector<Real> &x, Real &tol ) {
      const Real w0(1e5), h0(250);
      return -(horsepower(x[0],0)/h0 + warranty(x[1],0)/w0);
    }

    void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
      const Real w0(1e5), h0(250);
      g[0] = -horsepower(x[0],1)/h0;
      g[1] = -warranty(x[1],1)/w0;
    }
#if USE_HESSVEC
    void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real w0(1e5), h0(250);
      hv[0] = -horsepower(x[0],2)/h0 * v[0];
      hv[1] = -warranty(x[1],2)/w0 * v[1];
    }
#endif
  };

  template<class Real>
  class Constraint_CylinderHead : public StdConstraint<Real> {
  private:
    Real warranty(const Real flatness, const int deriv = 0) const {
      const Real w0(1e5), w1(1.5e4), w2(4), zero(0);
      return (deriv==0 ? w0 + w1*(w2 - flatness)
            :(deriv==1 ? -w1 : zero));
    }

    Real tcycle(const Real flatness, const int deriv = 0) const {
      const Real t0(45), t1(4.5), t2(4), tpwr(1.5), one(1);
      return (deriv==0 ? t0 + t1*std::pow(t2 - flatness, tpwr)
            :(deriv==1 ? -t1*tpwr*std::sqrt(t2 - flatness)
            : t1*tpwr*(tpwr-one)/std::sqrt(t2 - flatness)));
    }

    Real twall(const Real dintake, const int deriv = 0) const {
      const Real ointake(3.25), oexhaust(1.34), dexhaust(1.556), half(0.5), zero(0);
      return (deriv==0 ? ointake - oexhaust - half*(dintake + dexhaust)
            :(deriv==1 ? -half : zero));
    }

    Real smax(const Real tw, const int deriv = 0) const {
      const Real s0(750), spwr(2.5), one(1), two(2);
      return (deriv==0 ? s0 + one/std::pow(tw, spwr)
            :(deriv==1 ? -spwr/std::pow(tw, spwr+one)
            : spwr*(spwr+one)/std::pow(tw, spwr+two)));
    }

  public:
    Constraint_CylinderHead() {}

    void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
      const Real one(1), two(2), sixty(60), syield(3000), w0(1e5);
      c[0] = two*smax(twall(x[0],0),0)/syield - one;
      c[1] = one - warranty(x[1],0)/w0;
      c[2] = tcycle(x[1],0)/sixty - one;
    }

    void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real two(2), sixty(60), syield(3000), w0(1e5);
      jv[0] = two*smax(twall(x[0],0),1)/syield * twall(x[0],1) * v[0];
      jv[1] = -warranty(x[1],1)/w0 * v[1];
      jv[2] = tcycle(x[1],1)/sixty * v[1];
    }

    void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real two(2), sixty(60), syield(3000), w0(1e5);
      ajv[0] = two*smax(twall(x[0],0),1)/syield * twall(x[0],1) * v[0];
      ajv[1] = -warranty(x[1],1)/w0 * v[1] + tcycle(x[1],1)/sixty * v[2];
    }
#if USE_HESSVEC
    void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real two(2), sixty(60), syield(3000), w0(1e5);
      Real tw = twall(x[0],0);
      ahuv[0] = two*(smax(tw,2) * twall(x[0],1) * twall(x[0],1)
                + smax(tw,1) * twall(x[0],2))/syield * u[0] * v[0];
      ahuv[1] = (-warranty(x[1],2)/w0 * u[1] + tcycle(x[1],2)/sixty * u[2]) * v[1];
    }
#endif
  };

  template<class Real>
  class getCylinderHead : public TestProblem<Real> {
  public:
    getCylinderHead() {}

    Ptr<Objective<Real>> getObjective(void) const {
      return makePtr<Objective_CylinderHead<Real>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0)); 
      (*xp)[0] = static_cast<Real>(1.8);
      (*xp)[1] = static_cast<Real>(1.0);
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      (*xp)[0] = static_cast<Real>(2.122);
      (*xp)[1] = static_cast<Real>(1.769);
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> lp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      Ptr<std::vector<Real>> up    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      (*lp)[0] = static_cast<Real>(1.5);
      (*lp)[1] = static_cast<Real>(0.0);
      (*up)[0] = static_cast<Real>(2.164);
      (*up)[1] = static_cast<Real>(4.0);
      Ptr<Vector<Real>> l = makePtr<PrimalScaledStdVector<Real>>(lp,scale);
      Ptr<Vector<Real>> u = makePtr<PrimalScaledStdVector<Real>>(up,scale);
      return makePtr<Bounds<Real>>(l,u);
    }

    Ptr<Constraint<Real>> getInequalityConstraint(void) const {
      return makePtr<Constraint_CylinderHead<Real>>();
    }

    Ptr<Vector<Real>> getInequalityMultiplier(void) const {
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(3,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> lp    = makePtr<std::vector<Real>>(3,static_cast<Real>(0.0)); 
      return makePtr<DualScaledStdVector<Real>>(lp,scale);
    }

    Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(3,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> up    = makePtr<std::vector<Real>>(3,static_cast<Real>(0.0)); 
      Ptr<Vector<Real>> u = makePtr<DualScaledStdVector<Real>>(up,scale);
      return makePtr<Bounds<Real>>(*u,false);
    }
  };

}// End ZOO Namespace
}// End ROL Namespace

#endif
