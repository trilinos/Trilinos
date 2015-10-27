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
    \brief  Contains definitions for W. Hock and K. Schittkowski 25th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS25_HPP
#define ROL_HS25_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

  /** \brief W. Hock and K. Schittkowski 25th test function.
   */
  template<class Real>
  class Objective_HS25 : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;  
 
  typedef typename vector::size_type uint;

  private:
    std::vector<Real> u_vec_;
    uint u_size_;

  Teuchos::RCP<const vector> getVector( const V& x ) { 
    using Teuchos::dyn_cast;
    return dyn_cast<const SV>(x).getVector();  
  }

  Teuchos::RCP<vector> getVector( V& x ) {
    using Teuchos::dyn_cast;
    return dyn_cast<SV>(x).getVector();
  }

  public:
    Objective_HS25() {
      u_size_ = 99;
      for ( uint i = 0; i < u_size_; i++ ) {
        u_vec_.push_back(25.0 + std::pow((-50.0*std::log(0.01*(Real)(i+1))),2.0/3.0));
      }
    }

    Real value( const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);

      Real val = 0.0, f = 0.0;
      for ( uint i = 0; i < u_size_; i++ ) {
        f = -0.01*(Real)(i+1) + std::exp(-1.0/(*ex)[0] * std::pow(u_vec_[i]-(*ex)[1],(*ex)[2]));
        val += f*f;
      }
      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      RCP<vector> eg = getVector(g);
      g.zero();

      Real f = 0.0, df1 = 0.0, df2 = 0.0, df3 = 0.0, tmp = 0.0;
      for ( uint i = 0; i < u_size_; i++ ) {
        tmp = std::pow(u_vec_[i]-(*ex)[1],(*ex)[2])/(*ex)[0];
        f   = -0.01*(Real)(i+1) + std::exp(-tmp);
        df1 = std::exp(-tmp)*tmp/(*ex)[0];
        df2 = std::exp(-tmp)*(*ex)[2]*std::pow(u_vec_[i]-(*ex)[1],(*ex)[2]-1.0)/(*ex)[0];
        df3 = std::exp(-tmp)*tmp*std::log(u_vec_[i]-(*ex)[1]);
        (*eg)[0] += 2.0*f*df1;
        (*eg)[1] += 2.0*f*df2;
        (*eg)[2] += 2.0*f*df3;
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
 
      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      RCP<const vector> ev = getVector(v);
      RCP<vector> ehv = getVector(hv);
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x); 
      RCP<const vector> ev = getVector(v);
      RCP<vector> ehv = getVector(hv);
      hv.zero();
    }
  };

  template<class Real>
  void getHS25( Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<BoundConstraint<Real> > &con, 
                Vector<Real> &x0, Vector<Real> &x ) {
 
    typedef std::vector<Real>    vector;
    typedef ROL::Vector<Real>    V;
    typedef ROL::StdVector<Real> SV;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::dyn_cast;
    
    // Cast Initial Guess and Solution Vectors
    RCP<vector> x0p = dyn_cast<SV>(x0).getVector();
    RCP<vector> xp  = dyn_cast<SV>(x).getVector();

    uint n = xp->size();
    // Resize Vectors
    n = 3;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = rcp( new Objective_HS25<Real> );
    // Instantiate BoundConstraint

    RCP<vector> l_rcp = rcp( new vector(n,0.0) );
    RCP<vector> u_rcp = rcp( new vector(n,0.0) );
    l_rcp->push_back(0.1);
    l_rcp->push_back(0.0);
    l_rcp->push_back(0.0);
    u_rcp->push_back(100.0);
    u_rcp->push_back(25.6);
    u_rcp->push_back(5.0);

    RCP<V> l = rcp( new SV(l_rcp) );
    RCP<V> u = rcp( new SV(u_rcp) );
 
    con = rcp( new BoundConstraint<Real>(l,u) );

    // Get Initial Guess
    (*x0p)[0] = 100.0;
    (*x0p)[1] = 12.5;
    (*x0p)[2] = 3.0;

    // Get Solution
    (*xp)[0] = 50.0;
    (*xp)[1] = 25.0;
    (*xp)[2] = 1.5;
  }


} // End ZOO Namespace
} // End ROL Namespace

#endif
