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
 *  \brief Contains definitions for W. Hock and K. Schittkowski 32nd test problem
 *         which contains only inequality constraints.
 */

#ifndef ROL_HS29_HPP
#define ROL_HS29_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"


namespace ROL {
namespace ZOO {

template<class Real> 
class Objective_HS29 : public Objective<Real> {

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

    return -(*xp)[0]*(*xp)[1]*(*xp)[2];

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    
    Ptr<const vector> xp = getVector(x);
    Ptr<vector> gp = getVector(g);

    (*gp)[0] = -(*xp)[1]*(*xp)[2];
    (*gp)[1] = -(*xp)[0]*(*xp)[2];
    (*gp)[2] = -(*xp)[0]*(*xp)[1];

  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    
    Ptr<const vector> xp = getVector(x);
    Ptr<const vector> vp = getVector(v);
    Ptr<vector> hvp = getVector(hv);

    (*hvp)[0] = -( (*xp)[2]*(*vp)[1] + (*xp)[1]*(*vp)[2] );
    (*hvp)[1] = -( (*xp)[2]*(*vp)[0] + (*xp)[0]*(*vp)[2] );
    (*hvp)[2] = -( (*xp)[1]*(*vp)[0] + (*xp)[0]*(*vp)[1] );
 
  }

}; // class Objective_HS29


template<class Real>
class InequalityConstraint_HS29 : public Constraint<Real> {
  
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

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    

    Ptr<vector> cp = getVector(c);
    Ptr<const vector> xp = getVector(x); 

    (*cp)[0] = -std::pow((*xp)[0],2) - 2*std::pow((*xp)[1],2) - 4*std::pow((*xp)[2],2) + 48;

  }
 
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, 
                      const Vector<Real> &x, Real &tol ) {

    
    
    Ptr<vector> jvp = getVector(jv);
    Ptr<const vector> vp = getVector(v);
    Ptr<const vector> xp = getVector(x); 

    (*jvp)[0] = -2*(*xp)[0]*(*vp)[0] - 4*(*xp)[1]*(*vp)[1] - 8*(*xp)[2]*(*vp)[2];
       
  }
   
  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    

    Ptr<vector> ajvp = getVector(ajv);
    Ptr<const vector> vp = getVector(v);
    Ptr<const vector> xp = getVector(x);

    (*ajvp)[0] = -2*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] = -4*(*xp)[1]*(*vp)[0];
    (*ajvp)[2] = -8*(*xp)[2]*(*vp)[0];
 
  }

  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u,
                            const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    
  
    Ptr<vector> ahuvp = getVector(ahuv);
    Ptr<const vector> up = getVector(u);
    Ptr<const vector> vp = getVector(v);
    Ptr<const vector> xp = getVector(x);

    (*ahuvp)[0] = -2*(*up)[0]*(*vp)[0];
    (*ahuvp)[1] = -4*(*up)[0]*(*vp)[1];
    (*ahuvp)[2] = -8*(*up)[0]*(*vp)[2];
 
  }

}; // class InequalityConstraint_HS29


template<class Real>
class getHS29 : public TestProblem<Real> {
public:
  getHS29(void) {}

  Ptr<Objective<Real> > getObjective( void ) const {
    return makePtr<Objective_HS29<Real>>();
  }

  Ptr<Constraint<Real> > getInequalityConstraint( void ) const {
    return makePtr<InequalityConstraint_HS29<Real>>();
  }

  Ptr<BoundConstraint<Real> > getBoundConstraint( void ) const {
    // No Lower bound
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(3, ROL_NINF<Real>());
    
    // No upper bound
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(3, ROL_INF<Real>());
   
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
  
    return makePtr<Bounds<Real>>(l,u);
  }

  Ptr<Vector<Real> > getInitialGuess( void ) const {
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(3);
    (*x0p)[0] = 1.0;
    (*x0p)[1] = 1.0;
    (*x0p)[2] = 1.0;
  
    return makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real> > getSolution( void ) const {
    Ptr<std::vector<Real> > xp = makePtr<std::vector<Real>>(3);
    (*xp)[0] = 4.0;
    (*xp)[1] = 2.0*std::sqrt(2.0);
    (*xp)[2] = 2.0;
  
    return makePtr<StdVector<Real>>(xp);
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


#endif // ROL_HS29_HPP
