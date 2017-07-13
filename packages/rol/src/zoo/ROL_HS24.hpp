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
 *  \brief Contains definitions for W. Hock and K. Schittkowski 24th test problem
 *         which contains bound and inequality constraints.
 */

#ifndef ROL_HS24_HPP
#define ROL_HS24_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_InequalityConstraint.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {
namespace ZOO {

template<class Real>
class Objective_HS24 : public Objective<Real> {

  typedef std::vector<Real>   vector;
  typedef Vector<Real>        V;
  typedef StdVector<Real>     SV;

private:
  const Real rt3_;

public:

  Objective_HS24() : rt3_(std::sqrt(3)) {}

  Real value( const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector(); 

    return rt3_*(*xp)[0]*std::pow((*xp)[1],3)*((*xp)[0]-6)/81.0;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector(); 

    Teuchos::RCP<vector> gp = Teuchos::dyn_cast<SV>(g).getVector(); 

    (*gp)[0] = 2*rt3_*std::pow((*xp)[1],3)*((*xp)[0]-3)/81.0;
    (*gp)[1] = rt3_*(*xp)[0]*std::pow((*xp)[1],2)*((*xp)[0]-6)/27.0;

  }  

  
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector(); 
    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector(); 
    Teuchos::RCP<vector> hvp = Teuchos::dyn_cast<SV>(hv).getVector(); 

    Real a00 = pow((*xp)[1],3)/81.0;
    Real a01 = pow((*xp)[1],2)*((*xp)[0]-3)/27.0;
    Real a11 = (*xp)[1]*(std::pow((*xp)[0]-3,2)-9)/27.0;

    (*hvp)[0] = a00*(*vp)[0] + a01*(*vp)[1];
    (*hvp)[1] = a01*(*vp)[0] + a11*(*vp)[1]; 
    hv.scale(2*rt3_);

  }
}; // class Objective_HS24


template<class Real> 
class InequalityConstraint_HS24 : public InequalityConstraint<Real>  {

  typedef std::vector<Real>   vector;
  typedef Vector<Real>        V;
  typedef StdVector<Real>     SV;

private:
  const Real rt3_;

public:
  InequalityConstraint_HS24() : rt3_(std::sqrt(3)) {}

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> xp = Teuchos::dyn_cast<const SV>(x).getVector();
    Teuchos::RCP<vector> cp = Teuchos::dyn_cast<SV>(c).getVector();

    (*cp)[0] =  (*xp)[0]/rt3_ -      (*xp)[1];
    (*cp)[1] =  (*xp)[0]      + rt3_*(*xp)[1];
    (*cp)[2] = -(*xp)[0]      - rt3_*(*xp)[1] + 6;

  }

  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                      const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<vector> jvp = Teuchos::dyn_cast<SV>(jv).getVector();

    (*jvp)[0] =  (*vp)[0]/rt3_ -      (*vp)[1];
    (*jvp)[1] =  (*vp)[0]      + rt3_*(*vp)[1];
    (*jvp)[2] = -(*vp)[0]      - rt3_*(*vp)[1];


  }

  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, 
                             const Vector<Real> &x, Real &tol ) {

    Teuchos::RCP<const vector> vp = Teuchos::dyn_cast<const SV>(v).getVector();
    Teuchos::RCP<vector> ajvp = Teuchos::dyn_cast<SV>(ajv).getVector();

    (*ajvp)[0] = rt3_*(*vp)[0]/3 + (*vp)[1] - (*vp)[2];
    (*ajvp)[1] = -(*vp)[0] + rt3_*(*vp)[1] - rt3_*(*vp)[2];

  }   

  
  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u,
                            const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ahuv.zero();
  }

}; // class InequalityConstraint_HS32


template<class Real> 
Teuchos::RCP<Objective<Real> > getObjective_HS24( void ) {
  return Teuchos::rcp( new Objective_HS24<Real>() );
}

template<class Real> 
Teuchos::RCP<InequalityConstraint<Real> > getInequalityConstraint_HS24( void ) {
  return Teuchos::rcp( new InequalityConstraint_HS24<Real>() );
}


template<class Real> 
Teuchos::RCP<BoundConstraint<Real> > getBoundConstraint_HS24( void ) {

  // Lower bound is zero  
  Teuchos::RCP<std::vector<Real> > lp = Teuchos::rcp( new std::vector<Real>(2,0.0) );
  
  // No upper bound
  Teuchos::RCP<std::vector<Real> > up = Teuchos::rcp( new std::vector<Real>(2,ROL_INF<Real>()) );
 
  Teuchos::RCP<Vector<Real> > l = Teuchos::rcp( new StdVector<Real>(lp) );
  Teuchos::RCP<Vector<Real> > u = Teuchos::rcp( new StdVector<Real>(up) );

  return Teuchos::rcp( new BoundConstraint<Real>(l,u) );

}

template<class Real> 
Teuchos::RCP<Vector<Real> > getInitialGuess_HS24( void ) {

  Teuchos::RCP<std::vector<Real> > x0p = Teuchos::rcp( new std::vector<Real>(2) );
  (*x0p)[0] = 1.0;
  (*x0p)[1] = 0.5;

  return Teuchos::rcp( new StdVector<Real>(x0p) );
}

template<class Real> 
Teuchos::RCP<Vector<Real> > getSolution_HS24( void ) {

  Teuchos::RCP<std::vector<Real> > xp = Teuchos::rcp( new std::vector<Real>(2) );
  (*xp)[0] = 3.0;
  (*xp)[1] = std::sqrt(3.0);

  return Teuchos::rcp( new StdVector<Real>(xp) );
}

template<class Real> 
Teuchos::RCP<Vector<Real> > getInequalityMultiplier_HS24( void ) {
  
  Teuchos::RCP<std::vector<Real> > lp = Teuchos::rcp( new std::vector<Real>(3,0.0) );
  return Teuchos::rcp( new StdVector<Real>(lp) );

}




} // namespace ZOO
} // namespace ROL

#endif // ROL_HS24_HPP

