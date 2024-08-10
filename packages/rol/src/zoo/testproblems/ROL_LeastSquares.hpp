// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Least squares function.
 */
template<class Real>
class Objective_LeastSquares : public Objective<Real> {
  typedef typename std::vector<Real>::size_type uint; 

public:

  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    uint n    = ex->size();
    Real h   = 1.0/((Real)n+1.0);
    Real val = 0.0;
    Real res = 0.0;
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        res = 2.0*h*(5.0/6.0) + 1.0/h*((*ex)[i+1]-2.0*(*ex)[i]);
      }
      else if ( i == n-1 ) {
        res = 2.0*h*(5.0/6.0) + 1.0/h*((*ex)[i-1]-2.0*(*ex)[i]);
      }
      else {
        res = 2.0*h + 1.0/h*((*ex)[i-1]-2.0*(*ex)[i]+(*ex)[i+1]);
      }
      val += 0.5*res*res;
    }
    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > eg
      = dynamic_cast<StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    uint n  = ex->size();
    Real h = 1.0/((Real)n+1.0);
    std::vector<Real> res(n,0.0);
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*ex)[i+1]-2.0*(*ex)[i]);
      }
      else if ( i == n-1 ) {
        res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*ex)[i-1]-2.0*(*ex)[i]);
      }
      else {
        res[i] = 2.0*h + 1.0/h*((*ex)[i-1]-2.0*(*ex)[i]+(*ex)[i+1]);
      }
    }

    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        (*eg)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
      }
      else if ( i == n-1 ) {
        (*eg)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
      }
      else {
        (*eg)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
      }
    }
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > ehv
      = dynamic_cast<StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > ev
      = dynamic_cast<const StdVector<Real>&>(v).getVector();
    ROL::Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    uint n  = ex->size();
    Real h = 1.0/((Real)n+1.0);
    std::vector<Real> res(n,0.0);
    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        res[i] = 1.0/h*((*ev)[i+1]-2.0*(*ev)[i]);
      }
      else if ( i == n-1 ) {
        res[i] = 1.0/h*((*ev)[i-1]-2.0*(*ev)[i]);
      }
      else {
        res[i] = 1.0/h*((*ev)[i-1]-2.0*(*ev)[i]+(*ev)[i+1]);
      }
    }

    for (uint i=0; i<n; i++) {
      if ( i == 0 ) {
        (*ehv)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
      }
      else if ( i == n-1 ) {
        (*ehv)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
      }
      else {
        (*ehv)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
      }
    }
  }
#endif
};

template<class Real>
class getLeastSquares : public TestProblem<Real> {
public:
  getLeastSquares(void){}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_LeastSquares<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 32;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) {
      (*x0p)[i] = 0.0;
    }
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 32;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    Real h = 1.0/((Real)n+1.0), pt = 0.0;
    for( int i = 0; i < n; i++ ) {
      pt = (Real)(i+1)*h;
      (*xp)[i] = pt*(1.0-pt);
    }
    return ROL::makePtr<StdVector<Real>>(xp);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
