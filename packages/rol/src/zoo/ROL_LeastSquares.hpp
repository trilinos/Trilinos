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
    \brief  Contains definitions for least squares function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_LEASTSQUARES_HPP
#define ROL_LEASTSQUARES_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Least squares function.
   */
  template<class Real>
  class Objective_LeastSquares : public Objective<Real> {

    typedef std::vector<Real>  vector;
    typedef Vector<Real>       V;
    typedef StdVector<Real>    SV;     

    typedef typename vector::size_type uint; 

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

    Real value( const Vector<Real> &x, Real &tol ) {
      using Teuchos::RCP;
      RCP<const vector> xp = getVector(x);

      uint n    = xp->size();
      Real h   = 1.0/((Real)n+1.0);
      Real val = 0.0;
      Real res = 0.0;
      for (uint i=0; i<n; i++) {
        if ( i == 0 ) {
          res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
        }
        else if ( i == n-1 ) {
          res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
        }
        else {
          res = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);
        }
        val += 0.5*res*res;
      }

      return val;
   }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> xp = getVector(x);
      RCP<vector> gp = getVector(g);

      uint n  = xp->size();
      Real h = 1.0/((Real)n+1.0);
      vector res(n,0.0);
      for (uint i=0; i<n; i++) {
        if ( i == 0 ) {
          res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
        }
        else if ( i == n-1 ) {
          res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
        }
        else {
          res[i] = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);
        }
      }

      for (uint i=0; i<n; i++) {
        if ( i == 0 ) {
          (*gp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
        }
        else if ( i == n-1 ) {
          (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
        }
        else {
          (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
        }
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
 
      using Teuchos::RCP;
      RCP<const vector> xp = getVector(x);
      RCP<const vector> vp = getVector(v);
      RCP<vector> hvp = getVector(hv);

      uint n  = xp->size();
      Real h = 1.0/((Real)n+1.0);
      vector res(n,0.0);
      for (uint i=0; i<n; i++) {
        if ( i == 0 ) {
          res[i] = 1.0/h*((*vp)[i+1]-2.0*(*vp)[i]);
        }
        else if ( i == n-1 ) {
          res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]);
        }
        else {
          res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]+(*vp)[i+1]);
        }
      }

      for (uint i=0; i<n; i++) {
        if ( i == 0 ) {
          (*hvp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
        }
        else if ( i == n-1 ) {
          (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
        }
        else {
          (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
        }
      }
    }
#endif
  };

  template<class Real>
  void getLeastSquares( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {

    typedef std::vector<Real> vector;
    typedef StdVector<Real>   SV;

    typedef typename vector::size_type uint;
 
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::dyn_cast;

    // Cast Initial Guess and Solution Vectors
    RCP<vector> x0p = dyn_cast<SV>(x0).getVector();
    RCP<vector> xp  = dyn_cast<SV>(x).getVector();

    uint n = xp->size();

    // Resize Vectors
    n = 32;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = rcp( new Objective_LeastSquares<Real> );
    // Get Initial Guess
    for (uint i=0; i<n; i++) {
      (*x0p)[i] = 0.0;
    }
    // Get Solution
    Real h  = 1.0/((Real)n+1.0);
    Real pt = 0.0;
    for( uint i=0; i<n; i++ ) {
      pt = (Real)(i+1)*h;
      (*xp)[i] = pt*(1.0-pt);
    }
  }

} // End ZOO Namespace
} // End ROL Namespace

#endif
