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
    \brief  Contains definitions for W. Hock and K. Schittkowski 38th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS38_HPP
#define ROL_HS38_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

  /** \brief W. Hock and K. Schittkowski 38th test function.
   */
  template<class Real>
  class Objective_HS38 : public Objective<Real> {
   
  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;
  typedef StdVector<Real>   SV;  

  private:
  
    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const SV&>(x).getVector();
    }

    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<SV&>(x).getVector();  
    }

  public:
    Objective_HS38(void) {}

    Real value( const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> ex = getVector(x);
      return 100.0 * std::pow((*ex)[1] - std::pow((*ex)[0],2.0),2.0) + std::pow(1.0-(*ex)[0],2.0) + 
              90.0 * std::pow((*ex)[3] - std::pow((*ex)[2],2.0),2.0) + std::pow(1.0-(*ex)[2],2.0) +
              10.1 * (std::pow((*ex)[1] - 1.0,2.0) + std::pow((*ex)[3]-1.0,2.0)) + 
              19.8 * ((*ex)[1] - 1.0) * ((*ex)[3] - 1.0);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
 
      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<vector> eg = getVector(g); 

      (*eg)[0] = -4.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) * (*ex)[0] - 2.0 * (1.0-(*ex)[0]);
      (*eg)[1] = 2.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) + 
                 2.0 * 10.1 * ((*ex)[1] - 1.0) + 19.8*((*ex)[3] - 1.0); 
      (*eg)[2] = -4.0 * 90.0 * ((*ex)[3] - std::pow((*ex)[2],2.0)) * (*ex)[2] - 2.0 * (1.0-(*ex)[2]);
      (*eg)[3] = 2.0 * 90.0 * ((*ex)[3] - std::pow((*ex)[2],2.0)) + 
                 2.0 * 10.1 * ((*ex)[3] - 1.0) + 19.8*((*ex)[1] - 1.0); 
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<const vector> ev = getVector(v);
      ROL::Ptr<vector> ehv = getVector(hv);

      Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
      Real h12 = -4.0 * 100.0 * (*ex)[0];
      Real h13 = 0.0;
      Real h14 = 0.0;
      Real h21 = -4.0 * 100.0 * (*ex)[0];
      Real h22 = 2.0 * 100.0 + 2.0 * 10.1;
      Real h23 = 0.0;
      Real h24 = 19.8;
      Real h31 = 0.0;
      Real h32 = 0.0;
      Real h33 = -4.0 * 90.0 * (*ex)[3] + 12.0 * 90.0 * std::pow((*ex)[2],2.0) + 2.0; 
      Real h34 = -4.0 * 90.0 * (*ex)[2];
      Real h41 = 0.0;
      Real h42 = 19.8;
      Real h43 = -4.0 * 90.0 * (*ex)[2];
      Real h44 = 2.0 * 90.0 + 2.0 * 10.1;

      (*ehv)[0] = h11 * (*ev)[0] + h12 * (*ev)[1] + h13 * (*ev)[2] + h14 * (*ev)[3];
      (*ehv)[1] = h21 * (*ev)[0] + h22 * (*ev)[1] + h23 * (*ev)[2] + h24 * (*ev)[3];
      (*ehv)[2] = h31 * (*ev)[0] + h32 * (*ev)[1] + h33 * (*ev)[2] + h34 * (*ev)[3];
      (*ehv)[3] = h41 * (*ev)[0] + h42 * (*ev)[1] + h43 * (*ev)[2] + h44 * (*ev)[3];
    } 
#endif
  };

template<class Real>
class getHS38 : public TestProblem<Real> {
public:
  getHS38(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_HS38<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 4;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*x0p)[0] = -3.0; (*x0p)[1] = -1.0;
    (*x0p)[2] = -3.0; (*x0p)[3] = -1.0;
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 4;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = 1.0; (*xp)[1] = 1.0;
    (*xp)[2] = 1.0; (*xp)[3] = 1.0;
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    // Problem dimension
    int n = 4;
    // Instantiate BoundConstraint
    ROL::Ptr<std::vector<Real> > lp = ROL::makePtr<std::vector<Real>>(n,-10.0);
    ROL::Ptr<Vector<Real> > l = ROL::makePtr<StdVector<Real>>(lp);
    ROL::Ptr<std::vector<Real> > up = ROL::makePtr<std::vector<Real>>(n, 10.0);
    ROL::Ptr<Vector<Real> > u = ROL::makePtr<StdVector<Real>>(up);
    return ROL::makePtr<Bounds<Real>>(l,u);
  }
};


} // End ZOO Namespace
} // End ROL Namespace

#endif
