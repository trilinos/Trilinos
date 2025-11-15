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

/*! \file  test_19.cpp
    \brief Test ReducedConstraintSimOpt class

*/

#include "ROL_StdVector.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_Reduced_Constraint_SimOpt.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"
#include <cassert>

template<typename Real>
class constraint1 : public ROL::Constraint_SimOpt<Real> {
public:
  constraint1() {}
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(c.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> cs = dynamic_cast<ROL::StdVector<Real>&>(c);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(cs.getVector()))[0] = std::exp(z1*u1)-z2*z2;
  }
  void solve(ROL::Vector<Real> &c, ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(c.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> us = dynamic_cast<ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(us.getVector()))[0] = static_cast<Real>(2)*std::log(std::abs(z2)) / z1;
    constraint1<Real>::value(c,u,z,tol);
  }
  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(jv.dimension()==1);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> jvs = dynamic_cast<ROL::StdVector<Real>&>(jv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(jvs.getVector()))[0] = z1*std::exp(z1*u1)*v1;
  }
  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(jv.dimension()==1);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> jvs = dynamic_cast<ROL::StdVector<Real>&>(jv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(jvs.getVector()))[0] = u1*std::exp(z1*u1)*v1 - static_cast<Real>(2)*z2*v2;
  }
  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ijv.dimension()==1);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ijvs = dynamic_cast<ROL::StdVector<Real>&>(ijv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(ijvs.getVector()))[0] = v1 / (z1*std::exp(z1*u1));
  }
  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    constraint1<Real>::applyJacobian_1(ajv,v,u,z,tol);
  }
  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ajv.dimension()==2);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ajvs = dynamic_cast<ROL::StdVector<Real>&>(ajv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(ajvs.getVector()))[0] = u1*std::exp(z1*u1)*v1;
    (*(ajvs.getVector()))[1] = -static_cast<Real>(2)*z2*v1;
  }
  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &iajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    constraint1<Real>::applyInverseJacobian_1(iajv,v,u,z,tol);
  }
  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==1);
    assert(w.dimension()==1);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(ahwvs.getVector()))[0] = z1*z1*std::exp(z1*u1)*v1*w1;
  }
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==2);
    assert(w.dimension()==1);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(ahwvs.getVector()))[0] = std::exp(z1*u1)*(static_cast<Real>(1)+u1*z1)*v1*w1;
    (*(ahwvs.getVector()))[1] = static_cast<Real>(0);
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==1);
    assert(w.dimension()==1);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(ahwvs.getVector()))[0] = std::exp(z1*u1)*(static_cast<Real>(1)+u1*z1)*v1*w1;
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==2);
    assert(w.dimension()==1);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    (*(ahwvs.getVector()))[0] = u1*u1*std::exp(z1*u1)*v1*w1;
    (*(ahwvs.getVector()))[1] = -static_cast<Real>(2)*v2*w1;
  }
};


template<typename Real>
class constraint2 : public ROL::Constraint_SimOpt<Real> {
public:
  constraint2() {}
  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(c.dimension()==3);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> cs = dynamic_cast<ROL::StdVector<Real>&>(c);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(cs.getVector()))[0] = z1*z2*u1;
    (*(cs.getVector()))[1] = (z1-z2)*u1;
    (*(cs.getVector()))[2] = u1*u1;
  }
  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(jv.dimension()==3);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> jvs = dynamic_cast<ROL::StdVector<Real>&>(jv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    const Real two(2);
    Real v1 = (*(vs.getVector()))[0];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(jvs.getVector()))[0] = z1*z2*v1;
    (*(jvs.getVector()))[1] = (z1-z2)*v1;
    (*(jvs.getVector()))[2] = two*u1*v1;
  }
  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(jv.dimension()==3);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> jvs = dynamic_cast<ROL::StdVector<Real>&>(jv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(jvs.getVector()))[0] = z2*u1*v1 + z1*u1*v2;
    (*(jvs.getVector()))[1] = (v1-v2)*u1;
    (*(jvs.getVector()))[2] = static_cast<Real>(0);
  }
  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ajv.dimension()==1);
    assert(v.dimension()==3);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ajvs = dynamic_cast<ROL::StdVector<Real>&>(ajv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    const Real two(2);
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real v3 = (*(vs.getVector()))[2];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(ajvs.getVector()))[0] = z1*z2*v1 + (z1-z2)*v2 + two*u1*v3;
  }
  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ajv.dimension()==2);
    assert(v.dimension()==3);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ajvs = dynamic_cast<ROL::StdVector<Real>&>(ajv);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real u1 = (*(us.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(ajvs.getVector()))[0] = (z2*u1*v1 + u1*v2);
    (*(ajvs.getVector()))[1] = (z1*u1*v1 - u1*v2);
  }
  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==1);
    assert(w.dimension()==3);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    const Real two(2);
    Real w3 = (*(ws.getVector()))[2];
    Real v1 = (*(vs.getVector()))[0];
    (*(ahwvs.getVector()))[0] = two*v1*w3;
  }
  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==2);
    assert(w.dimension()==3);
    assert(v.dimension()==1);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real w2 = (*(ws.getVector()))[1];
    Real v1 = (*(vs.getVector()))[0];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(ahwvs.getVector()))[0] = (z2*v1*w1 + v1*w2);
    (*(ahwvs.getVector()))[1] = (z1*v1*w1 - v1*w2);
  }
  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==1);
    assert(w.dimension()==3);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real w2 = (*(ws.getVector()))[1];
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real z1 = (*(zs.getVector()))[0];
    Real z2 = (*(zs.getVector()))[1];
    (*(ahwvs.getVector()))[0] = (v1*z2+z1*v2)*w1 + (v1-v2)*w2;
  }
  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    assert(ahwv.dimension()==2);
    assert(w.dimension()==3);
    assert(v.dimension()==2);
    assert(u.dimension()==1);
    assert(z.dimension()==2);
    ROL::StdVector<Real> ahwvs = dynamic_cast<ROL::StdVector<Real>&>(ahwv);
    const ROL::StdVector<Real> ws = dynamic_cast<const ROL::StdVector<Real>&>(w);
    const ROL::StdVector<Real> vs = dynamic_cast<const ROL::StdVector<Real>&>(v);
    const ROL::StdVector<Real> us = dynamic_cast<const ROL::StdVector<Real>&>(u);
    const ROL::StdVector<Real> zs = dynamic_cast<const ROL::StdVector<Real>&>(z);
    Real w1 = (*(ws.getVector()))[0];
    Real v1 = (*(vs.getVector()))[0];
    Real v2 = (*(vs.getVector()))[1];
    Real u1 = (*(us.getVector()))[0];
    (*(ahwvs.getVector()))[0] = v2*u1*w1;
    (*(ahwvs.getVector()))[1] = v1*u1*w1;
  }
};

int main(int argc, char *argv[]) {
  using RealT = double;

  ROL::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

//  RealT errtol = std::sqrt(ROL::ROL_THRESHOLD<RealT>());

  int errorFlag  = 0;

  // *** Test body.

  try {

    unsigned c1_dim = 1; // Constraint1 dimension
    unsigned c2_dim = 3; // Constraint1 dimension
    unsigned u_dim  = 1; // State dimension
    unsigned z_dim  = 2; // Control dimension

    auto c1  = ROL::makePtr<ROL::StdVector<RealT>>(c1_dim);
    auto c2  = ROL::makePtr<ROL::StdVector<RealT>>(c2_dim);
    auto u   = ROL::makePtr<ROL::StdVector<RealT>>(u_dim);
    auto z   = ROL::makePtr<ROL::StdVector<RealT>>(z_dim);
    auto vc1 = ROL::makePtr<ROL::StdVector<RealT>>(c1_dim);
    auto vc2 = ROL::makePtr<ROL::StdVector<RealT>>(c2_dim);
    auto vu  = ROL::makePtr<ROL::StdVector<RealT>>(u_dim);
    auto vz  = ROL::makePtr<ROL::StdVector<RealT>>(z_dim);
    auto du  = ROL::makePtr<ROL::StdVector<RealT>>(u_dim);
    auto dz  = ROL::makePtr<ROL::StdVector<RealT>>(z_dim);
    c1->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    c2->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    u->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    z->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    vc1->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    vc2->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    vu->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    vz->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    du->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    dz->randomize(static_cast<RealT>(-1),static_cast<RealT>(1));
    
    auto con1 = ROL::makePtr<constraint1<RealT>>();
    auto con2 = ROL::makePtr<constraint2<RealT>>();
    auto stateStore = ROL::makePtr<ROL::VectorController<RealT>>();
    auto rcon = ROL::makePtr<ROL::Reduced_Constraint_SimOpt<RealT>>(con2,con1,stateStore,u,z,vc1,c2,true,false);

    con1->checkSolve(*u,*z,*c1,true,*outStream);
    con1->checkAdjointConsistencyJacobian_1(*vc1,*vu,*u,*z,true,*outStream);
    con1->checkAdjointConsistencyJacobian_2(*vc1,*vz,*u,*z,true,*outStream);
    con1->checkInverseJacobian_1(*c1,*vu,*u,*z,true,*outStream);
    con1->checkInverseAdjointJacobian_1(*c1,*vu,*u,*z,true,*outStream);
    con1->checkApplyJacobian_1(*u,*z,*vu,*vc1,true,*outStream);
    con1->checkApplyJacobian_2(*u,*z,*vz,*vc1,true,*outStream);
    con1->checkApplyAdjointHessian_11(*u,*z,*vc1,*vu,*du,true,*outStream);
    con1->checkApplyAdjointHessian_12(*u,*z,*vc1,*vu,*dz,true,*outStream);
    con1->checkApplyAdjointHessian_21(*u,*z,*vc1,*vz,*du,true,*outStream);
    con1->checkApplyAdjointHessian_22(*u,*z,*vc1,*vz,*dz,true,*outStream);

    con2->checkAdjointConsistencyJacobian_1(*vc2,*vu,*u,*z,true,*outStream);
    con2->checkAdjointConsistencyJacobian_2(*vc2,*vz,*u,*z,true,*outStream);
    con2->checkApplyJacobian_1(*u,*z,*vu,*vc2,true,*outStream);
    con2->checkApplyJacobian_2(*u,*z,*vz,*vc2,true,*outStream);
    con2->checkApplyAdjointHessian_11(*u,*z,*vc2,*vu,*du,true,*outStream);
    con2->checkApplyAdjointHessian_12(*u,*z,*vc2,*vu,*dz,true,*outStream);
    con2->checkApplyAdjointHessian_21(*u,*z,*vc2,*vz,*du,true,*outStream);
    con2->checkApplyAdjointHessian_22(*u,*z,*vc2,*vz,*dz,true,*outStream);

    rcon->checkAdjointConsistencyJacobian(*vc2,*vz,*z,true,*outStream);
    rcon->checkApplyJacobian(*z,*vz,*vc2,true,*outStream);
    rcon->checkApplyAdjointHessian(*z,*vc2,*vz,*dz,true,*outStream);
  }   
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;


}

