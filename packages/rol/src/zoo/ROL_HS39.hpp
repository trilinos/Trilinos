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
    \brief  Contains definitions for W. Hock and K. Schittkowski 39th test function.
    \author Created by G. von Winckel
 */


#ifndef ROL_HS39_HPP
#define ROL_HS39_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint_Partitioned.hpp"
#include "ROL_BoundConstraint.hpp"
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
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();
    return -(*xp)[0];
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();
    Teuchos::RCP<vector> gp = Teuchos::dyn_cast<SV>(g).getVector();

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
class EqualityConstraint_HS39a : public EqualityConstraint<Real> {

  typedef std::vector<Real> vector;
  typedef StdVector<Real>   SV;

public:
 
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<vector> cp = Teuchos::dyn_cast<SV>(c).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*cp)[0] = (*xp)[1]-std::pow((*xp)[0],3)-std::pow((*xp)[2],2);
  }  

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {

    Teuchos::RCP<vector> jvp = Teuchos::dyn_cast<SV>(jv).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*jvp)[0] = (*vp)[1] - 3.0*(*xp)[0]*(*xp)[0]*(*vp)[0] - 2.0*(*xp)[2]*(*vp)[2];

  }

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<vector> ajvp = Teuchos::dyn_cast<SV>(ajv).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*ajvp)[0] = -3.0*(*xp)[0]*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] = (*vp)[0];
    (*ajvp)[2] = -2.0*(*xp)[2]*(*vp)[0];
    (*ajvp)[3] = 0.0;
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                           const Vector<Real> &v, const Vector<Real> &x,
                           Real &tol) {

    Teuchos::RCP<vector> ahuvp = Teuchos::dyn_cast<SV>(ahuv).getVector();
    Teuchos::RCP<const vector> up = Teuchos::dyn_cast<const SV>(u).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*ahuvp)[0] = -6.0*(*up)[0]*(*xp)[0]*(*vp)[0]; 
    (*ahuvp)[1] = 0.0;
    (*ahuvp)[2] = -2.0*(*up)[0]*(*vp)[2];
    (*ahuvp)[3] = 0.0;

  }


};

// Second of two equality constraints
template<class Real>
class EqualityConstraint_HS39b : public EqualityConstraint<Real> {

  typedef std::vector<Real> vector;
  typedef StdVector<Real>   SV;

public:
 
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
    Teuchos::RCP<vector> cp = Teuchos::dyn_cast<SV>(c).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*cp)[0] = std::pow((*xp)[0],2)-(*xp)[1]-std::pow((*xp)[3],2);
  }  

  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) {

    Teuchos::RCP<vector> jvp = Teuchos::dyn_cast<SV>(jv).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*jvp)[0] = 2.0*(*xp)[0]*(*vp)[0] - (*vp)[1] - 2.0*(*xp)[3]*(*vp)[3];

  }

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v,
                             const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<vector> ajvp = Teuchos::dyn_cast<SV>(ajv).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

    (*ajvp)[0] = 2.0*(*xp)[0]*(*vp)[0];
    (*ajvp)[1] = -(*vp)[0]; 
    (*ajvp)[2] = 0.0; 
    (*ajvp)[3] = -2.0*(*vp)[0]*(*xp)[3];
  }

  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                           const Vector<Real> &v, const Vector<Real> &x,
                           Real &tol) {

    Teuchos::RCP<vector> ahuvp = Teuchos::dyn_cast<SV>(ahuv).getVector();
    Teuchos::RCP<const vector> up = Teuchos::dyn_cast<const SV>(u).getVector();
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();

  //  (*cp)[0] = std::pow((*xp)[0],2)-(*xp)[1]-std::pow((*xp)[3],2);

    (*ahuvp)[0] = 2.0*(*up)[0]*(*vp)[0];
    (*ahuvp)[1] = 0.0;
    (*ahuvp)[2] = 0.0;
    (*ahuvp)[3] = -2.0*(*up)[0]*(*vp)[3];
  }


};




} // End ZOO Namespace
} // End ROL Namespace

#endif
