// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_10.cpp
    \brief Show how to use CompositeConstraint interface.
*/

#include "ROL_CompositeConstraint_SimOpt.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

template<class Real>
class valConstraint : public ROL::Constraint_SimOpt<Real> {
public:
  valConstraint(void) : ROL::Constraint_SimOpt<Real>() {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp
      = dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();

    Real half(0.5), two(2);
    // C(0) = U(0) - Z(0)
    (*cp)[0] = (*up)[0]-(*zp)[0];
    // C(1) = 0.5 * (U(0) + U(1) - Z(0))^2
    (*cp)[1] = half*std::pow((*up)[0]+(*up)[1]-(*zp)[0],two);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp
      = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*jvp)[0] = (*vp)[0];
    (*jvp)[1] = ((*up)[0] + (*up)[1] - (*zp)[0]) * ((*vp)[0] + (*vp)[1]);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp
      = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*jvp)[0] = -(*vp)[0];
    (*jvp)[1] = ((*zp)[0] - (*up)[0] - (*up)[1]) * (*vp)[0];
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ajvp
      = dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ajvp)[0] = (*vp)[0] + ((*up)[0] + (*up)[1] - (*zp)[0]) * (*vp)[1];
    (*ajvp)[1] = ((*up)[0] + (*up)[1] - (*zp)[0]) * (*vp)[1];
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ajvp
      = dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ajvp)[0] = ((*zp)[0] - (*up)[0] - (*up)[1]) * (*vp)[1] - (*vp)[0];
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = (*wp)[1] * ((*vp)[0] + (*vp)[1]);
    (*ahwvp)[1] = (*wp)[1] * ((*vp)[0] + (*vp)[1]);
  }

  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = -(*wp)[1] * ((*vp)[0] + (*vp)[1]);
  }

  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = -(*wp)[1] * (*vp)[0];
    (*ahwvp)[1] = -(*wp)[1] * (*vp)[0];
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = (*wp)[1] * (*vp)[0];
  }
};

template<class Real>
class redConstraint : public ROL::Constraint_SimOpt<Real> {
public:
  redConstraint(void) : ROL::Constraint_SimOpt<Real>() {}

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp
      = dynamic_cast<ROL::StdVector<Real>&>(c).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();

    const Real one(1), two(2);
    // C = exp(U) - (Z^2 + 1)
    (*cp)[0] = std::exp((*up)[0])-(std::pow((*zp)[0],two) + one);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp
      = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*jvp)[0] = std::exp((*up)[0]) * (*vp)[0];
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > jvp
      = dynamic_cast<ROL::StdVector<Real>&>(jv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();

    const Real two(2);
    (*jvp)[0] = -two * (*zp)[0] * (*vp)[0];
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ajvp
      = dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ajvp)[0] = std::exp((*up)[0]) * (*vp)[0];
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &ajv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ajvp
      = dynamic_cast<ROL::StdVector<Real>&>(ajv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();

    const Real two(2);
    (*ajvp)[0] = -two * (*zp)[0] * (*vp)[0];
  }

  void applyInverseJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ijvp
      = dynamic_cast<ROL::StdVector<Real>&>(ijv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ijvp)[0] = (*vp)[0] / std::exp((*up)[0]);
  }

  void applyInverseAdjointJacobian_1(ROL::Vector<Real> &ijv, const ROL::Vector<Real> &v,
                                     const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ijvp
      = dynamic_cast<ROL::StdVector<Real>&>(ijv).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ijvp)[0] = (*vp)[0] / std::exp((*up)[0]);
  }

  void applyAdjointHessian_11(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = std::exp((*up)[0]) * (*wp)[0] * (*vp)[0];
  }

  void applyAdjointHessian_12(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = static_cast<Real>(0);
  }

  void applyAdjointHessian_21(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = static_cast<Real>(0);
  }

  void applyAdjointHessian_22(ROL::Vector<Real> &ahwv, const ROL::Vector<Real> &w, const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > ahwvp
      = dynamic_cast<ROL::StdVector<Real>&>(ahwv).getVector();
    ROL::Ptr<const std::vector<Real> > wp
      = dynamic_cast<const ROL::StdVector<Real>&>(w).getVector();
    ROL::Ptr<const std::vector<Real> > vp
      = dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > up
      = dynamic_cast<const ROL::StdVector<Real>&>(u).getVector();
    ROL::Ptr<const std::vector<Real> > zp
      = dynamic_cast<const ROL::StdVector<Real>&>(z).getVector();
    (*ahwvp)[0] = static_cast<Real>(-2) * (*wp)[0] * (*vp)[0];
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {

    int dim = 2;
    int dimz = 1;
    ROL::Ptr<std::vector<RealT> > ustd  = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT> > dustd = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT> > zstd  = ROL::makePtr<std::vector<RealT>>(dimz);
    ROL::Ptr<std::vector<RealT> > dzstd = ROL::makePtr<std::vector<RealT>>(dimz);
    ROL::Ptr<std::vector<RealT> > cstd  = ROL::makePtr<std::vector<RealT>>(dim);
    ROL::Ptr<std::vector<RealT> > czstd = ROL::makePtr<std::vector<RealT>>(dimz);
    ROL::Ptr<std::vector<RealT> > sstd  = ROL::makePtr<std::vector<RealT>>(dimz);
    ROL::Ptr<std::vector<RealT> > dsstd = ROL::makePtr<std::vector<RealT>>(dimz);

    (*ustd)[0]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*ustd)[1]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*dustd)[0] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*dustd)[1] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*zstd)[0]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*dzstd)[0] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*cstd)[0]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*cstd)[1]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*czstd)[0] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*sstd)[0]  = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    (*dsstd)[0] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);

    ROL::Ptr<ROL::Vector<RealT> > u  = ROL::makePtr<ROL::StdVector<RealT>>(ustd);
    ROL::Ptr<ROL::Vector<RealT> > du = ROL::makePtr<ROL::StdVector<RealT>>(dustd);
    ROL::Ptr<ROL::Vector<RealT> > z  = ROL::makePtr<ROL::StdVector<RealT>>(zstd);
    ROL::Ptr<ROL::Vector<RealT> > dz = ROL::makePtr<ROL::StdVector<RealT>>(dzstd);
    ROL::Ptr<ROL::Vector<RealT> > c  = ROL::makePtr<ROL::StdVector<RealT>>(cstd);
    ROL::Ptr<ROL::Vector<RealT> > cz = ROL::makePtr<ROL::StdVector<RealT>>(czstd);
    ROL::Ptr<ROL::Vector<RealT> > s  = ROL::makePtr<ROL::StdVector<RealT>>(sstd);
    ROL::Ptr<ROL::Vector<RealT> > ds = ROL::makePtr<ROL::StdVector<RealT>>(dsstd);

    ROL::Vector_SimOpt<RealT> x(u,s);
    ROL::Vector_SimOpt<RealT> dx(du,ds);
    ROL::Vector_SimOpt<RealT> y(s,z);
    ROL::Vector_SimOpt<RealT> dy(ds,dz);
    ROL::Vector_SimOpt<RealT> w(u,z);
    ROL::Vector_SimOpt<RealT> dw(du,dz);

    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > valCon = ROL::makePtr<valConstraint<RealT>>();
    valCon->checkAdjointConsistencyJacobian_1(*c,*du,*u,*s,true,*outStream);
    valCon->checkAdjointConsistencyJacobian_2(*c,*dz,*u,*s,true,*outStream);
    valCon->checkApplyJacobian_1(*u,*s,*du,*c,true,*outStream);
    valCon->checkApplyJacobian_2(*u,*s,*ds,*c,true,*outStream);
    valCon->checkApplyJacobian(x,dx,*c,true,*outStream);
    valCon->checkApplyAdjointHessian_11(*u,*s,*c,*du,*u,true,*outStream);
    valCon->checkApplyAdjointHessian_12(*u,*s,*c,*du,*s,true,*outStream);
    valCon->checkApplyAdjointHessian_21(*u,*s,*c,*ds,*u,true,*outStream);
    valCon->checkApplyAdjointHessian_22(*u,*s,*c,*ds,*s,true,*outStream);
    valCon->checkApplyAdjointHessian(x,*c,dx,x,true,*outStream);

    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > redCon = ROL::makePtr<redConstraint<RealT>>();
    redCon->checkAdjointConsistencyJacobian_1(*cz,*ds,*s,*z,true,*outStream);
    redCon->checkAdjointConsistencyJacobian_2(*cz,*dz,*s,*z,true,*outStream);
    redCon->checkInverseJacobian_1(*cz,*ds,*s,*z,true,*outStream); 
    redCon->checkInverseAdjointJacobian_1(*ds,*cz,*s,*z,true,*outStream); 
    redCon->checkApplyJacobian_1(*s,*z,*ds,*cz,true,*outStream);
    redCon->checkApplyJacobian_2(*s,*z,*dz,*cz,true,*outStream);
    redCon->checkApplyJacobian(y,dy,*cz,true,*outStream);
    redCon->checkApplyAdjointHessian_11(*s,*z,*cz,*ds,*s,true,*outStream);
    redCon->checkApplyAdjointHessian_12(*s,*z,*cz,*ds,*z,true,*outStream);
    redCon->checkApplyAdjointHessian_21(*s,*z,*cz,*dz,*s,true,*outStream);
    redCon->checkApplyAdjointHessian_22(*s,*z,*cz,*dz,*z,true,*outStream);
    redCon->checkApplyAdjointHessian(y,*cz,dy,y,true,*outStream);

    ROL::CompositeConstraint_SimOpt<RealT> con(valCon,redCon,*c,*cz,*u,*s,*z);
    con.checkAdjointConsistencyJacobian_1(*c,*du,*u,*z,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(*c,*dz,*u,*z,true,*outStream);
    con.checkApplyJacobian_1(*u,*z,*du,*c,true,*outStream);
    con.checkApplyJacobian_2(*u,*z,*dz,*c,true,*outStream);
    con.checkApplyJacobian(w,dw,*c,true,*outStream);
    con.checkApplyAdjointHessian_11(*u,*z,*c,*du,*u,true,*outStream);
    con.checkApplyAdjointHessian_12(*u,*z,*c,*du,*z,true,*outStream);
    con.checkApplyAdjointHessian_21(*u,*z,*c,*dz,*u,true,*outStream);
    con.checkApplyAdjointHessian_22(*u,*z,*c,*dz,*z,true,*outStream);
    con.checkApplyAdjointHessian(w,*c,dw,w,true,*outStream);
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

