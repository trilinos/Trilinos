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
 *         which contains both inequality and equality constraints.
 */

#ifndef ROL_HS32_HPP
#define ROL_HS32_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"


namespace ROL {
namespace ZOO {

template<class Real> 
class Objective_HS32 : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;

private:

  ROL::SharedPointer<const vector> getVector( const V& x ) { 
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::SharedPointer<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  
  Real value( const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<const vector> xp = getVector(x);

    Real term1 = (*xp)[0]+3*(*xp)[1]+(*xp)[2];
    Real term2 = (*xp)[0]-(*xp)[1];
    return term1*term1 + 4*term2*term2;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    
    
    ROL::SharedPointer<vector> gp = getVector(g);
    ROL::SharedPointer<const vector> xp = getVector(x);

    (*gp)[0] = 10*(*xp)[0] -  2*(*xp)[1] + 2*(*xp)[2];
    (*gp)[1] = -2*(*xp)[0] + 26*(*xp)[1] + 6*(*xp)[2];
    (*gp)[2] =  2*(*xp)[0] +  6*(*xp)[1] + 2*(*xp)[2];
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    
 
    ROL::SharedPointer<vector> hvp = getVector(hv);
    ROL::SharedPointer<const vector> vp = getVector(v);

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

  ROL::SharedPointer<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector(); 
  }

  ROL::SharedPointer<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  EqualityConstraint_HS32() {}

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> cp = getVector(c);
    ROL::SharedPointer<const vector> xp = getVector(x);

    (*cp)[0] = 1.0 - (*xp)[0] - (*xp)[1] - (*xp)[2];
  }

  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                      const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> jvp = getVector(jv);
    ROL::SharedPointer<const vector> vp = getVector(v);

    (*jvp)[0] = - (*vp)[0] - (*vp)[1] - (*vp)[2];

  } 

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    

    ROL::SharedPointer<vector> ajvp = getVector(ajv);
    ROL::SharedPointer<const vector> vp = getVector(v);
     
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
  InequalityConstraint_HS32(void) {}
  
  ROL::SharedPointer<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::SharedPointer<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> cp = getVector(c);
    ROL::SharedPointer<const vector> xp = getVector(x);

    (*cp)[0] = 6*(*xp)[1]+4*(*xp)[2]-std::pow((*xp)[0],3)-3.0;
    (*cp)[1] = (*xp)[0];
    (*cp)[2] = (*xp)[1];
    (*cp)[3] = (*xp)[2];

  }
 
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, 
                      const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> jvp = getVector(jv);
    ROL::SharedPointer<const vector> vp = getVector(v);
    ROL::SharedPointer<const vector> xp = getVector(x);

    (*jvp)[0] = -3*(*xp)[0]*(*xp)[0]*(*vp)[0]+6*(*vp)[1]+4*(*vp)[2];
    (*jvp)[1] = (*vp)[0];
    (*jvp)[2] = (*vp)[1];
    (*jvp)[3] = (*vp)[2];
        
  }
   
  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> ajvp = getVector(ajv);
    ROL::SharedPointer<const vector> vp = getVector(v); 
    ROL::SharedPointer<const vector> xp = getVector(x); 

    (*ajvp)[0] = -3*(*xp)[0]*(*xp)[0]*(*vp)[0] + (*vp)[1];
    (*ajvp)[1] =  6*(*vp)[0] + (*vp)[2]; 
    (*ajvp)[2] =  4*(*vp)[0] + (*vp)[3];
    
  }

  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u,
                            const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    
    ROL::SharedPointer<vector> ahuvp = getVector(ahuv); 
    ROL::SharedPointer<const vector> up = getVector(u);
    ROL::SharedPointer<const vector> vp = getVector(v);
    ROL::SharedPointer<const vector> xp = getVector(x);
 
    (*ahuvp)[0] = -6*(*up)[0]*(*vp)[0]*(*xp)[0];
    (*ahuvp)[1] = 0.0;
    (*ahuvp)[2] = 0.0;

  }

}; // class InequalityConstraint_HS32

}
} // namespace ROL


#endif // ROL_HS32_HPP
