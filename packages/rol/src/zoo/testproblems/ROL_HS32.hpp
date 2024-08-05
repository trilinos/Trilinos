// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
 *  \brief Contains definitions for W. Hock and K. Schittkowski 32nd test problem
 *         which contains both inequality and equality constraints.
 */

#ifndef ROL_HS32_HPP
#define ROL_HS32_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

template<class Real> 
class Objective_HS32 : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;

private:

  Ptr<const vector> getVector( const V& x ) { 
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  
  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const vector> xp = getVector(x);

    Real term1 = (*xp)[0]+3*(*xp)[1]+(*xp)[2];
    Real term2 = (*xp)[0]-(*xp)[1];
    return term1*term1 + 4*term2*term2;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<vector> gp = getVector(g);
    Ptr<const vector> xp = getVector(x);

    (*gp)[0] = 10*(*xp)[0] -  2*(*xp)[1] + 2*(*xp)[2];
    (*gp)[1] = -2*(*xp)[0] + 26*(*xp)[1] + 6*(*xp)[2];
    (*gp)[2] =  2*(*xp)[0] +  6*(*xp)[1] + 2*(*xp)[2];
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<vector> hvp = getVector(hv);
    Ptr<const vector> vp = getVector(v);

    (*hvp)[0] = 10*(*vp)[0] -  2*(*vp)[1] + 2*(*vp)[2];
    (*hvp)[1] = -2*(*vp)[0] + 26*(*vp)[1] + 6*(*vp)[2];
    (*hvp)[2] =  2*(*vp)[0] +  6*(*vp)[1] + 2*(*vp)[2];
  }

}; // class Objective_HS32


template<class Real>
class EqualityConstraint_HS32 : public Constraint<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;

private:

  Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector(); 
  }

  Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  EqualityConstraint_HS32() {}

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    const Real one(1);
    Ptr<vector> cp = getVector(c);
    Ptr<const vector> xp = getVector(x);

    (*cp)[0] = one - (*xp)[0] - (*xp)[1] - (*xp)[2];
  }

  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                      const Vector<Real> &x, Real &tol ) {
    Ptr<vector> jvp = getVector(jv);
    Ptr<const vector> vp = getVector(v);

    (*jvp)[0] = - (*vp)[0] - (*vp)[1] - (*vp)[2];
  } 

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {
    Ptr<vector> ajvp = getVector(ajv);
    Ptr<const vector> vp = getVector(v);
     
    (*ajvp)[0] = -(*vp)[0];
    (*ajvp)[1] = -(*vp)[0];
    (*ajvp)[2] = -(*vp)[0];
  }

  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u,
                            const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ahuv.zero();
  }

}; // class EqualityConstraint_HS32


template<class Real>
class InequalityConstraint_HS32 : public Constraint<Real> {

  typedef std::vector<Real> vector;  
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;

private:
  Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const SV&>(x).getVector();
  }

  Ptr<vector> getVector( V& x ) {
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  InequalityConstraint_HS32(void) {}

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    const Real three(3), four(4), six(6);
    Ptr<vector> cp = getVector(c);
    Ptr<const vector> xp = getVector(x);

    (*cp)[0] = six*(*xp)[1]+four*(*xp)[2]-std::pow((*xp)[0],three)-three;
  }
 
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, 
                      const Vector<Real> &x, Real &tol ) {
    const Real three(3), four(4), six(6);
    Ptr<vector> jvp = getVector(jv);
    Ptr<const vector> vp = getVector(v);
    Ptr<const vector> xp = getVector(x);

    (*jvp)[0] = -three*(*xp)[0]*(*xp)[0]*(*vp)[0]+six*(*vp)[1]+four*(*vp)[2];
  }
   
  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {
    const Real three(3), four(4), six(6);
    Ptr<vector> ajvp = getVector(ajv);
    Ptr<const vector> vp = getVector(v); 
    Ptr<const vector> xp = getVector(x); 

    (*ajvp)[0] = -three*(*xp)[0]*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] =  six*(*vp)[0];
    (*ajvp)[2] =  four*(*vp)[0];
  }

  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u,
                            const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    const Real zero(0), six(6);
    Ptr<vector> ahuvp = getVector(ahuv); 
    Ptr<const vector> up = getVector(u);
    Ptr<const vector> vp = getVector(v);
    Ptr<const vector> xp = getVector(x);
 
    (*ahuvp)[0] = -six*(*up)[0]*(*vp)[0]*(*xp)[0];
    (*ahuvp)[1] = zero;
    (*ahuvp)[2] = zero;
  }

}; // class Constraint_HS32


template<class Real>
class getHS32 : public TestProblem<Real> {
public:
  getHS32(void) {}

  Ptr<Objective<Real> > getObjective( void ) const {
    return makePtr<Objective_HS32<Real>>();
  }

  Ptr<Constraint<Real> > getEqualityConstraint( void ) const {
    return makePtr<EqualityConstraint_HS32<Real>>();
  }

  Ptr<Constraint<Real> > getInequalityConstraint( void ) const {
    return makePtr<InequalityConstraint_HS32<Real>>();
  }

  Ptr<BoundConstraint<Real> > getBoundConstraint( void ) const {
    // Lower bound zero
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(3, static_cast<Real>(0));
    
    // No upper bound
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(3, ROL_INF<Real>());
   
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
  
    return makePtr<Bounds<Real>>(l,u);
  }

  Ptr<Vector<Real> > getInitialGuess( void ) const {
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(3);
    (*x0p)[0] = 0.1;
    (*x0p)[1] = 0.7;
    (*x0p)[2] = 0.2;

    return makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real> > getSolution( const int i = 0 ) const {
    Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(3);
    (*xp)[0] = 0.0;
    (*xp)[1] = 0.0;
    (*xp)[2] = 1.0;
  
    return makePtr<StdVector<Real>>(xp);
  }

  Ptr<Vector<Real> > getEqualityMultiplier( void ) const {
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(1,0.0);
    return makePtr<StdVector<Real>>(lp);
  }

  Ptr<Vector<Real> > getInequalityMultiplier( void ) const {
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(1,0.0);
    return makePtr<StdVector<Real>>(lp);
  }

  Ptr<BoundConstraint<Real>> getSlackBoundConstraint(void) const {
    // Lower bound is zero  
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(1,0.0);
    
    // No upper bound
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(1,ROL_INF<Real>());
   
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
  
    return makePtr<Bounds<Real>>(l,u);
  }
};

}
} // namespace ROL


#endif // ROL_HS32_HPP
