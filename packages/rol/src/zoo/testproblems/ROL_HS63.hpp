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
    \brief  Contains definitions for W. Hock and K. Schittkowski 63rd test function.
    \author Created by D. P. Kouri
 */


#ifndef ROL_HS63_HPP
#define ROL_HS63_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Types.hpp"
#include "ROL_StdObjective.hpp"
#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 63rd test function.
 */

template<class Real>
class Objective_HS63 : public StdObjective<Real> {
public:
  Real value( const std::vector<Real> &x, Real &tol ) {
    const Real c1(1000), c2(2);
    return c1 - x[0]*x[0] - c2*x[1]*x[1] - x[2]*x[2] - x[0]*x[1] - x[0]*x[2];
  }

  void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    g[0] = -c2*x[0]-x[1]-x[2];
    g[1] = -c2*c2*x[1]-x[0];
    g[2] = -c2*x[2]-x[0];
  }

  void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    const Real c2(2);
    hv[0] = -c2*v[0]-v[1]-v[2];
    hv[1] = -c2*c2*v[1]-v[0];
    hv[2] = -c2*v[2]-v[0];
  } 
};

template<class Real>
class Constraint_HS63a : public StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(8), c2(14), c3(7), c4(56);
    c[0] = c1*x[0] + c2*x[1] + c3*x[2] - c4;
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c1(8), c2(14), c3(7);
    jv[0] = c1*v[0] + c2*v[1] + c3*v[2];
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c1(8), c2(14), c3(7);
    ajv[0] = c1*v[0];
    ajv[1] = c2*v[0];
    ajv[2] = c3*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }
};

template<class Real>
class Constraint_HS63b : public StdConstraint<Real> {
public:
  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) {
    const Real c1(25);
    c[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - c1;
  }

  void applyJacobian(std::vector<Real> &jv, const std::vector<Real> &v,
                     const std::vector<Real> &x, Real &tol) {
    const Real c1(2);
    jv[0] = c1*(x[0]*v[0] + x[1]*v[1] + x[2]*v[2]);
  }

  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v,
                             const std::vector<Real> &x, Real &tol ) {
    const Real c1(2);
    ajv[0] = c1*x[0]*v[0];
    ajv[1] = c1*x[1]*v[0];
    ajv[2] = c1*x[2]*v[0];
  }

  void applyAdjointHessian(std::vector<Real> &ahuv, const std::vector<Real> &u,
                           const std::vector<Real> &v, const std::vector<Real> &x,
                           Real &tol) {
    const Real c1(2);
    ahuv[0] = c1*v[0]*u[0];
    ahuv[1] = c1*v[1]*u[0];
    ahuv[2] = c1*v[2]*u[0];
  }
};

template<class Real>
class getHS63 : public TestProblem<Real> {
public:
  getHS63(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    return ROL::makePtr<Objective_HS63<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    int n = 3;
    return ROL::makePtr<StdVector<Real>>(n,static_cast<Real>(2));
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    int n = 3;
    ROL::Ptr<std::vector<Real>> xp = ROL::makePtr<std::vector<Real>>(n);
    (*xp)[0] = static_cast<Real>(3.512118414);
    (*xp)[1] = static_cast<Real>(0.2169881741);
    (*xp)[2] = static_cast<Real>(3.552174034);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    std::vector<Ptr<Constraint<Real>>> cvec(2);
    cvec[0] = makePtr<Constraint_HS63a<Real>>();
    cvec[1] = makePtr<Constraint_HS63b<Real>>();
    return ROL::makePtr<Constraint_Partitioned<Real>>(cvec);
  }

  Ptr<Vector<Real>> getEqualityMultiplier(void) const {
    std::vector<Ptr<Vector<Real>>> lvec(2);
    lvec[0] = makePtr<StdVector<Real>>(1,static_cast<Real>(0));
    lvec[1] = makePtr<StdVector<Real>>(1,static_cast<Real>(0));
    return ROL::makePtr<PartitionedVector<Real>>(lvec);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    int n = 3;
    Ptr<Vector<Real>> l = makePtr<StdVector<Real>>(n,static_cast<Real>(0.0));
    return makePtr<Bounds<Real>>(*l,true);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
