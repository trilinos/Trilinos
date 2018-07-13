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
    \brief  Contains definitions for W. Hock and K. Schittkowski 5th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS5_HPP
#define ROL_HS5_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

  /** \brief W. Hock and K. Schittkowski 5th test function.
   */
  template<class Real>
  class Objective_HS5 : public Objective<Real> {
  
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
    Objective_HS5(void) {}

    Real value( const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> ex = getVector(x);

      return std::sin((*ex)[0] + (*ex)[1]) + std::pow((*ex)[0]-(*ex)[1],2.0) - 1.5*(*ex)[0] + 2.5*(*ex)[1] + 1.0;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<vector> eg = getVector(g);

      (*eg)[0] = std::cos((*ex)[0] + (*ex)[1]) + 2.0*((*ex)[0]-(*ex)[1]) - 1.5;
      (*eg)[1] = std::cos((*ex)[0] + (*ex)[1]) - 2.0*((*ex)[0]-(*ex)[1]) + 2.5;;
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
 
      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<const vector> ev = getVector(v);
      ROL::Ptr<vector> ehv = getVector(hv);
  
      Real h11 = -std::sin((*ex)[0] + (*ex)[1]) + 2.0;
      Real h22 = -std::sin((*ex)[0] + (*ex)[1]) + 2.0;
      Real h12 = -std::sin((*ex)[0] + (*ex)[1]) - 2.0;
      Real h21 = -std::sin((*ex)[0] + (*ex)[1]) - 2.0;

      (*ehv)[0] = h11 * (*ev)[0] + h12 * (*ev)[1];
      (*ehv)[1] = h21 * (*ev)[0] + h22 * (*ev)[1];
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<const vector> ev = getVector(v);
      ROL::Ptr<vector> ehv = getVector(hv);
  
      Real h11 = -std::sin((*ex)[0] + (*ex)[1]) + 2.0;
      Real h22 = -std::sin((*ex)[0] + (*ex)[1]) + 2.0;
      Real h12 = -std::sin((*ex)[0] + (*ex)[1]) - 2.0;
      Real h21 = -std::sin((*ex)[0] + (*ex)[1]) - 2.0;

      (*ehv)[0] = 1.0/(h11*h22 - h12*h21) * (h22 * (*ev)[0] - h12 * (*ev)[1]);
      (*ehv)[1] = 1.0/(h11*h22 - h12*h21) * (-h21 * (*ev)[0] + h11 * (*ev)[1]);
    }
  };

template<class Real>
class getHS5 : public TestProblem<Real> {
public:
  getHS5(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_HS5<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 2;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*x0p)[0] = 0.0; (*x0p)[1] = 0.0;
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 2;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = -ROL::ScalarTraits<Real>::pi()/3.0 + 1.0/2.0; (*xp)[1] = -ROL::ScalarTraits<Real>::pi()/3.0 - 1.0/2.0;
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    // Problem dimension
    int n = 2;
    // Instantiate BoundConstraint
    ROL::Ptr<std::vector<Real> > lp = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*lp)[0] = -1.5; (*lp)[1] = -3.0;
    ROL::Ptr<Vector<Real> > l = ROL::makePtr<StdVector<Real>>(lp);
    ROL::Ptr<std::vector<Real> > up = ROL::makePtr<std::vector<Real>>(n,0.0);
    (*up)[0] = 4.0; (*up)[1] = 3.0;
    ROL::Ptr<Vector<Real> > u = ROL::makePtr<StdVector<Real>>(up);
    return ROL::makePtr<Bounds<Real>>(l,u);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
