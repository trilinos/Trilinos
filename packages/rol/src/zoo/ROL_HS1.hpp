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
    \brief  Contains definitions for W. Hock and K. Schittkowski 1st test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS1_HPP
#define ROL_HS1_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

  /** \brief W. Hock and K. Schittkowski 1st test function.
   */
  template<class Real>
  class Objective_HS1 : public Objective<Real> {

    typedef std::vector<Real> vector;
    typedef Vector<Real>      V;
    typedef StdVector<Real>   SV;  

  private:

    Teuchos::RCP<const vector> getVector( const V& x ) {
      using Teuchos::dyn_cast;
      return dyn_cast<const SV>(x).getVector();
    }

    Teuchos::RCP<vector> getVector( V& x ) {
      using Teuchos::dyn_cast;
      return dyn_cast<SV>(x).getVector();
    }    

  public:
    Objective_HS1(void) {}

    Real value( const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      return 100.0 * std::pow((*ex)[1] - std::pow((*ex)[0],2.0),2.0) + std::pow(1.0-(*ex)[0],2.0);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
     
      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      RCP<vector> eg = getVector(g);
      (*eg)[0] = -4.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) * (*ex)[0] - 2.0 * (1.0-(*ex)[0]);
      (*eg)[1] = 2.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)); 
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

     using Teuchos::RCP; 
     RCP<const vector> ex = getVector(x);
     RCP<const vector> ev = getVector(v);
     RCP<vector> ehv = getVector(hv);

      Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
      Real h22 = 2.0 * 100.0;
      Real h12 = -4.0 * 100.0 * (*ex)[0];
      Real h21 = -4.0 * 100.0 * (*ex)[0];

      (*ehv)[0] = h11 * (*ev)[0] + h12 * (*ev)[1];
      (*ehv)[1] = h21 * (*ev)[0] + h22 * (*ev)[1];
    } 
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      RCP<const vector> ev = getVector(v);
      RCP<vector> ehv = getVector(hv); 
      
      Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
      Real h22 = 2.0 * 100.0;
      Real h12 = -4.0 * 100.0 * (*ex)[0];
      Real h21 = -4.0 * 100.0 * (*ex)[0];
  
      (*ehv)[0] = 1.0/(h11*h22 - h12*h21) * (h22 * (*ev)[0] - h12 * (*ev)[1]);
      (*ehv)[1] = 1.0/(h11*h22 - h12*h21) * (-h21 * (*ev)[0] + h11 * (*ev)[1]);
    }
  };

  template<class Real>
  void getHS1( Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<BoundConstraint<Real> > &con, 
                Vector<Real> &x0, Vector<Real> &x ) {
 
    typedef std::vector<Real> vector;
    typedef Vector<Real>      V;
    typedef StdVector<Real>   SV;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::dyn_cast;
   
   // Cast Initial Guess and Solution Vectors
   RCP<vector> x0p = dyn_cast<SV>(x0).getVector();
   RCP<vector> xp  = dyn_cast<SV>(x).getVector();

   int n = xp->size();
   // Resize Vectors
   n = 2;
   x0p->resize(n);
   xp->resize(n);

   // Instantiate Objective Function
   obj = Teuchos::rcp( new Objective_HS1<Real> );

   // Instantiate BoundConstraint
   RCP<vector> l_rcp = rcp( new vector(n,0.0) ); 

   (*l_rcp)[0] = -ROL_OVERFLOW; 
   (*l_rcp)[1] = -1.5;

   RCP<vector> u_rcp = rcp( new vector(n,0.0) ); 
 
   (*u_rcp)[0] = ROL_OVERFLOW;  
   (*u_rcp)[1] = ROL_OVERFLOW;
    
   RCP<V> l = rcp( new SV(l_rcp) );
   RCP<V> u = rcp( new SV(u_rcp) );

   con = rcp( new BoundConstraint<Real>(l,u) );

   // Get Initial Guess
   (*x0p)[0] =  -2.0;
   (*x0p)[1] =  1.0;
   // Get Solution
   (*xp)[0] = 1.0;
   (*xp)[1] = 1.0;
 }


} // End ZOO Namespace
} // End ROL Namespace

#endif
