// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for W. Hock and K. Schittkowski 45th test function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_HS45_HPP
#define ROL_HS45_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

  /** \brief W. Hock and K. Schittkowski 45th test function.
   */
  template<class Real>
  class Objective_HS45 : public Objective<Real> {

    typedef std::vector<Real> vector;
    typedef Vector<Real>      V;
    typedef StdVector<Real>   SV;

    typedef typename vector::size_type uint;

  private: 
    uint  dim_;
    Real fact_;

    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const SV&>(x).getVector();  
    }

    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<SV&>(x).getVector();
    }

  public:
    Objective_HS45(int dim = 5) : dim_((uint)dim) {
      fact_ = 1.0;
      for ( uint i = 0; i < dim_; i++ ) {
        fact_ *= (Real)(i+1);
      } 
    }

    Real value( const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> ex = getVector(x); 
      Real prod = 1.0;
      for ( uint i = 0; i < dim_; i++ ) {
        prod *= (*ex)[i];
      }
      return 2.0 - prod/fact_;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<vector> eg = getVector(g);

      Real prod = 1.0;
      for ( uint j = 0; j < dim_; j++ ) {
        for ( uint i = 0; i < dim_; i++ ) {
          if ( j != i ) {
            prod *= (*ex)[i];
          }
        }
        (*eg)[j] = -prod/fact_;
        prod = 1.0;
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> ex = getVector(x);
      ROL::Ptr<const vector> ev = getVector(v);
      ROL::Ptr<vector> ehv = getVector(hv);

      hv.zero();
      Real prod = 1.0;
      for ( uint l = 0; l < dim_; l++ ) {
        for ( uint j = 0; j < dim_; j++ ) {
          if ( l != j ) {
            for ( uint i = 0; i < dim_; i++ ) {
              if ( j != i && l != i ) { 
                prod *= (*ex)[i];
              }
            }
            (*ehv)[l] += -prod/fact_*(*ev)[j];
          }
          prod = 1.0;
        }
      }
    } 
#endif
  };

template<class Real>
class getHS45 : public TestProblem<Real> {
public:
  getHS45(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Problem dimension
    int n = 10;
    // Instantiate Objective Function
    return ROL::makePtr<Objective_HS45<Real>>(n);
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 10;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) {
      (*x0p)[i] = 2.0;
    }
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 10;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) {
      (*xp)[i] = (Real)(i+1);
    }
    return ROL::makePtr<StdVector<Real>>(xp);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    // Problem dimension
    int n = 10;
    // Instantiate BoundConstraint
    ROL::Ptr<std::vector<Real> > lp = ROL::makePtr<std::vector<Real>>(n,0.0);
    ROL::Ptr<std::vector<Real> > up = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) { 
      (*lp)[i] = 0.0;
      (*up)[i] = static_cast<Real>(i+1);
    }
    ROL::Ptr<Vector<Real> > l = ROL::makePtr<StdVector<Real>>(lp);
    ROL::Ptr<Vector<Real> > u = ROL::makePtr<StdVector<Real>>(up);
    return ROL::makePtr<Bounds<Real>>(l,u);
  }
};


} // End ZOO Namespace
} // End ROL Namespace

#endif
