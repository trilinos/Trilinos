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
    \brief  Contains definitions for the discrete boundary value problem (More, Garbow, Hillstrom, 1981).
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_BVP_HPP
#define ROL_BVP_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {
namespace ZOO {

  /** \brief The discrete boundary value problem.
   */
  template<class Real>
  class Objective_BVP : public Objective<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;  

  typedef typename vector::size_type uint;

  private: 
    uint dim_;

    Teuchos::RCP<const vector> getVector( const V& x ) {
      using Teuchos::dyn_cast;
      using Teuchos::getConst;
      return dyn_cast<const SV>(getConst(x)).getVector();
    }

    Teuchos::RCP<vector> getVector( V& x ) {
      using Teuchos::dyn_cast;
      return dyn_cast<SV>(x).getVector();
    }


  public:
    Objective_BVP(void) : dim_(20) {}

    Real value( const Vector<Real> &x, Real &tol ) {
 
      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      Real val = 0.0;
      Real f   = 0.0;
      Real h   = 1.0/((Real)(dim_) + 1.0);
      for ( uint i = 0; i < dim_; i++ ) {
        f = 2.0*(*ex)[i] + h*h*std::pow((*ex)[i] + (Real)(i+1)*h + 1.0,3.0)/2.0; 
        if ( i < (dim_-1) ) { f -= (*ex)[i+1]; } 
        if ( i > 0 )              { f -= (*ex)[i-1]; }
        val += f*f;
      }
      return val; 
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;    
      RCP<const vector> ex = getVector(x);
      RCP<vector> eg = getVector(g);
      g.zero();
      Real h  = 1.0/((Real)(dim_) + 1.0);
      vector f(dim_,0.0);

      for ( uint i = 0; i < dim_; i++ ) {
        f[i] = 2.0*(*ex)[i] + h*h*std::pow((*ex)[i] + (Real)(i+1)*h + 1.0,3.0)/2.0;
        if ( i < (dim_-1) ) { f[i] -= (*ex)[i+1]; }
        if ( i > 0)               { f[i] -= (*ex)[i-1]; }
      }
      Real df = 0.0;
      for ( uint i = 0; i < dim_; i++ ) {
        df = (2.0 + 3.0*h*h*std::pow((*ex)[i] + (Real)(i+1)*h + 1.0,2.0)/2.0)*f[i];
        if ( i < (dim_-1) ) { df -= f[i+1]; }
        if ( i > 0 )        { df -= f[i-1]; }
        (*eg)[i] += 2.0*df;
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      using Teuchos::RCP;
      RCP<const vector> ex = getVector(x);
      RCP<const vector> ev = getVector(v);
      RCP<vector> ehv = getVector(hv);

      hv.zero();
      Real h = 1.0/((Real)(dim_) + 1.0);
      Real f = 0.0, df = 0.0, dfn = 0.0, hf = 0.0;
      for ( uint i = 0; i < dim_; i++ ) {
        f  = 2.0*(*ex)[i] + h*h*std::pow((*ex)[i] + (Real)(i+1)*h + 1.0,3.0)/2.0;
        df = 2.0 + 3.0/2.0 * h*h * std::pow((*ex)[i] + (Real)(i+1)*h + 1.0,2.0);
        hf = 3.0 * h*h * ((*ex)[i] + (Real)(i+1)*h + 1.0);
        if ( i < (dim_-2) ) {
          (*ehv)[i] += 2.0*(*ev)[i+2];
        }
        if ( i < (dim_-1) ) {
          f -= (*ex)[i+1];
          dfn = 2.0 + 3.0/2.0 * h*h * std::pow((*ex)[i+1] + (Real)(i+2)*h + 1.0,2.0);
          (*ehv)[i] -= 2.0*(df + dfn)*(*ev)[i+1];
          (*ehv)[i] += 2.0*(*ev)[i];
        }
        if ( i > 0 ) {
          f -= (*ex)[i-1];
          dfn = 2.0 + 3.0/2.0 * h*h * std::pow((*ex)[i-1] + (Real)(i)*h + 1.0,2.0);
          (*ehv)[i] -= 2.0*(df + dfn)*(*ev)[i-1];
          (*ehv)[i] += 2.0*(*ev)[i];
        }  
        if ( i > 1 ) {
          (*ehv)[i] += 2.0*(*ev)[i-2];
        } 
        (*ehv)[i] += 2.0*(hf*f + df*df)*(*ev)[i];  
      }
    } 
#endif
  };

  template<class Real>
  void getBVP( Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<BoundConstraint<Real> > &con, 
                Vector<Real> &x0, Vector<Real> &x ) {

    typedef std::vector<Real>     vector;
    typedef Vector<Real>          V;
    typedef StdVector<Real>       SV; 

    typedef typename vector::size_type uint;
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::dyn_cast;

    // Cast Initial Guess and Solution Vectors
    RCP<vector> x0p = dyn_cast<SV>(x0).getVector(); 
    RCP<vector> xp  = dyn_cast<SV>(x).getVector();
 
    uint n = xp->size();
    // Resize Vectors
    n = 20;
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_BVP<Real> );

    // Instantiate BoundConstraint

    RCP<vector> l_rcp = rcp( new vector() );
    RCP<vector> u_rcp = rcp( new vector() );

    RCP<V> l = rcp( new SV(l_rcp) );
    RCP<V> u = rcp( new SV(u_rcp) );

    std::vector<Real> val(n,0.0); 
    val[0] = 0.1*0.2321;
    val[1] = -0.1*0.4520;
    val[2] = -0.1*0.6588;
    val[3] = -0.1*0.8514;
    val[4] = -0.1*1.0288;
    val[5] = -0.1*1.1985;
    val[6] = -0.1*1.3322;
    val[7] = -0.1*1.4553;
    val[8] = -0.1*1.5571;
    val[9] = -0.1*1.6354;
    val[10] = -0.1*1.6881;
    val[11] = -0.1*1.7127;
    val[12] = -0.1*1.7060;
    val[13] = -0.1*1.6650;
    val[14] = -0.1*1.5856;
    val[15] = -0.1*1.4636;
    val[16] = -0.1*1.2938;
    val[17] = -0.1*1.0702;
    val[18] = -0.1*0.7858;
    val[19] = -0.1*0.4323;
    for ( uint i = 0; i < n; i++ ) { 
      if ( i%2 == 0 ) {  
        l_rcp->push_back(std::max(-0.2*(Real)(n),val[i]+0.1));
        u_rcp->push_back(std::min( 0.2*(Real)(n),val[i]+1.1));
      }
      else {
        l_rcp->push_back(-0.2*(Real)(n));
        u_rcp->push_back( 0.2*(Real)(n));
      }
    }
    
    con = rcp( new BoundConstraint<Real>(l,u) );
    // Get Initial Guess
    Real h = 1.0/((Real)n + 1.0);
    for ( uint i = 0; i < n; i++ ) {
      (*x0p)[i] = (Real)(i+1)*h*((Real)(i+1)*h - 1.0);
    }
    con->project(x0);
    // Get Solution
    (*xp)[0]  = 1.2321000000000001e-01; 
    (*xp)[1]  = 2.1743122909175336e-01;
    (*xp)[2]  = 2.8625218549543746e-01;
    (*xp)[3]  = 3.3309751851140840e-01;
    (*xp)[4]  = 3.6117201714254760e-01;
    (*xp)[5]  = 3.7342787212179440e-01;
    (*xp)[6]  = 3.7255212003706123e-01;
    (*xp)[7]  = 3.6096984201471016e-01;
    (*xp)[8]  = 3.4085861052124522e-01;
    (*xp)[9]  = 3.1417024791439530e-01;
    (*xp)[10] = 2.8265678244892922e-01; 
    (*xp)[11] = 2.4789833165179542e-01;
    (*xp)[12] = 2.1133139591375166e-01;
    (*xp)[13] = 1.7427666644258599e-01;
    (*xp)[14] = 1.3796594229036069e-01;
    (*xp)[15] = 1.0356813245768780e-01;
    (*xp)[16] = 7.2214621084083663e-02;
    (*xp)[17] = 4.5024529114833199e-02;
    (*xp)[18] = 2.3130648161534966e-02;
    (*xp)[19] = 7.7070870882527927e-03;
  }


}// End ZOO Namespace
}// End ROL Namespace

#endif
