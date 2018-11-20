// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for a quartic test problem.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_QUARTIC_HPP
#define ROL_QUARTIC_HPP

#include "ROL_ScaledStdVector.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

  template<class Real>
  class Objective_Quartic : public StdObjective<Real> {
  public:
    Objective_Quartic() {}

    Real value( const std::vector<Real> &x, Real &tol ) {
      const Real one(1);
      return std::pow(x[0]-one, 4) + std::pow(x[1]-one, 4);
    }

    void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
      const Real one(1), four(4);
      g[0] = four * std::pow(x[0]-one, 3);
      g[1] = four * std::pow(x[1]-one, 3);
    }
#if USE_HESSVEC
    void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real one(1), twelve(12);
      hv[0] = twelve * std::pow(x[0]-one, 2) * v[0];
      hv[1] = twelve * std::pow(x[1]-one, 2) * v[1];
    }
#endif
  };

  template<class Real>
  class Constraint_Quartic : public StdConstraint<Real> {
  public:
    Constraint_Quartic() {}

    void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
      const Real half(0.5);
      c[0] = std::pow(x[0], 2) - half*x[1];
      c[1] = std::pow(x[1], 2) - half*x[0];
    }

    void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real half(0.5), two(2);
      jv[0] = two*x[0]*v[0] - half*v[1];
      jv[1] = two*x[1]*v[1] - half*v[0];
    }

    void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real half(0.5), two(2);
      ajv[0] = two*x[0]*v[0] - half*v[1];
      ajv[1] = two*x[1]*v[1] - half*v[0];
    }
#if USE_HESSVEC
    void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
      const Real two(2);
      ahuv[0] = two*u[0]*v[0];
      ahuv[1] = two*u[1]*v[1];
    }
#endif
  };

  template<class Real>
  class getQuartic : public TestProblem<Real> {
  public:
    getQuartic() {}

    Ptr<Objective<Real>> getObjective(void) const {
      return makePtr<Objective_Quartic<Real>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0)); 
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> xp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.5));
      return makePtr<PrimalScaledStdVector<Real>>(xp,scale);
    }

    Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
      int n = 2;
      Ptr<std::vector<Real>> scale = makePtr<std::vector<Real>>(n,static_cast<Real>(1.0));
      Ptr<std::vector<Real>> lp    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      Ptr<std::vector<Real>> up    = makePtr<std::vector<Real>>(n,static_cast<Real>(0.0));
      (*lp)[0] = static_cast<Real>( 0.5);
      (*lp)[1] = static_cast<Real>(-2.9);
      (*up)[0] = static_cast<Real>( 5.8);
      (*up)[1] = static_cast<Real>( 2.9);
      Ptr<Vector<Real>> l = makePtr<PrimalScaledStdVector<Real>>(lp,scale);
      Ptr<Vector<Real>> u = makePtr<PrimalScaledStdVector<Real>>(up,scale);
      return makePtr<Bounds<Real>>(l,u);
    }

    Ptr<Constraint<Real>> getInequalityConstraint(void) const {
      return makePtr<Constraint_Quartic<Real>>();
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
