// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 39th test function.
    \author Created by G. von Winckel
 */


#ifndef ROL_HS39_HPP
#define ROL_HS39_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 39th test function.
 *
 *  Exact solution x* = (1, 0, 0, 0)
 *  f(x*) = -1
 */

template<class Real>
class Objective_HS39 : public Objective<Real> {

  typedef std::vector<Real>    vector;

  typedef StdVector<Real>   SV;



public:

  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();
    return -(*xp)[0];
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();
    ROL::Ptr<vector> gp = dynamic_cast<SV&>(g).getVector();

    (*gp)[0] = -1.0;
    (*gp)[1] =  0.0;
    (*gp)[2] =  0.0;
    (*gp)[3] =  0.0; 

  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    hv.zero();  
  } 
};

// First of two equality constraints
template<class Real>
class Constraint_HS39a : public Constraint<Real> {

  typedef std::vector<Real> vector;
  typedef StdVector<Real>   SV;

public:
  Constraint_HS39a(void) {}
 
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<vector> cp = dynamic_cast<SV&>(c).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*cp)[0] = (*xp)[1]-std::pow((*xp)[0],3)-std::pow((*xp)[2],2);
  }  

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {

    ROL::Ptr<vector> jvp = dynamic_cast<SV&>(jv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*jvp)[0] = (*vp)[1] - 3.0*(*xp)[0]*(*xp)[0]*(*vp)[0] - 2.0*(*xp)[2]*(*vp)[2];

  }

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<vector> ajvp = dynamic_cast<SV&>(ajv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*ajvp)[0] = -3.0*(*xp)[0]*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] = (*vp)[0];
    (*ajvp)[2] = -2.0*(*xp)[2]*(*vp)[0];
    (*ajvp)[3] = 0.0;
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                           const Vector<Real> &v, const Vector<Real> &x,
                           Real &tol) {

    ROL::Ptr<vector> ahuvp = dynamic_cast<SV&>(ahuv).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*ahuvp)[0] = -6.0*(*up)[0]*(*xp)[0]*(*vp)[0]; 
    (*ahuvp)[1] = 0.0;
    (*ahuvp)[2] = -2.0*(*up)[0]*(*vp)[2];
    (*ahuvp)[3] = 0.0;

  }


};

// Second of two equality constraints
template<class Real>
class Constraint_HS39b : public Constraint<Real> {

  typedef std::vector<Real> vector;
  typedef StdVector<Real>   SV;

public:
  Constraint_HS39b(void) {}
 
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<vector> cp = dynamic_cast<SV&>(c).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*cp)[0] = std::pow((*xp)[0],2)-(*xp)[1]-std::pow((*xp)[3],2);
  }  

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {

    ROL::Ptr<vector> jvp = dynamic_cast<SV&>(jv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*jvp)[0] = 2.0*(*xp)[0]*(*vp)[0] - (*vp)[1] - 2.0*(*xp)[3]*(*vp)[3];

  }

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    ROL::Ptr<vector> ajvp = dynamic_cast<SV&>(ajv).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

    (*ajvp)[0] = 2.0*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] = -(*vp)[0]; 
    (*ajvp)[2] = 0.0; 
    (*ajvp)[3] = -2.0*(*vp)[0]*(*xp)[3];
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                           const Vector<Real> &v, const Vector<Real> &x,
                           Real &tol) {

    ROL::Ptr<vector> ahuvp = dynamic_cast<SV&>(ahuv).getVector();
    ROL::Ptr<const vector> up = dynamic_cast<const SV&>(u).getVector();
    ROL::Ptr<const vector> vp = dynamic_cast<const SV&>(v).getVector();
    ROL::Ptr<const vector> xp = dynamic_cast<const SV&>(x).getVector();

  //  (*cp)[0] = std::pow((*xp)[0],2)-(*xp)[1]-std::pow((*xp)[3],2);

    (*ahuvp)[0] = 2.0*(*up)[0]*(*vp)[0];
    (*ahuvp)[1] = 0.0;
    (*ahuvp)[2] = 0.0;
    (*ahuvp)[3] = -2.0*(*up)[0]*(*vp)[3];
  }


};

template<class Real>
class getHS39 : public TestProblem<Real> {
public:
  getHS39(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_HS39<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 4;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,2.0);
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 4;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = static_cast<Real>(1);
    (*xp)[1] = static_cast<Real>(1);
    (*xp)[2] = static_cast<Real>(0);
    (*xp)[3] = static_cast<Real>(0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<Constraint<Real>> getEqualityConstraint(void) const {
    std::vector<Ptr<Constraint<Real>>> cvec(2);
    cvec[0] = makePtr<Constraint_HS39a<Real>>();
    cvec[1] = makePtr<Constraint_HS39b<Real>>();
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
