// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "ROL_TestProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace ZOO {

/** \brief W. Hock and K. Schittkowski 1st test function.
 */
template<class Real>
class Objective_HS1 : public Objective<Real> {
public:
  Objective_HS1(void) {}

  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();
    return 100.0 * std::pow((*ex)[1] - std::pow((*ex)[0],2.0),2.0) + std::pow(1.0-(*ex)[0],2.0);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > eg
      = dynamic_cast<StdVector<Real>&>(g).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();
   
    (*eg)[0] = -4.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)) * (*ex)[0] - 2.0 * (1.0-(*ex)[0]);
    (*eg)[1] = 2.0 * 100.0 * ((*ex)[1] - std::pow((*ex)[0],2.0)); 
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > ehv
      = dynamic_cast<StdVector<Real>&>(hv).getVector();
    Ptr<const std::vector<Real> > ev
      = dynamic_cast<const StdVector<Real>&>(v).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();

    Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
    Real h22 = 2.0 * 100.0;
    Real h12 = -4.0 * 100.0 * (*ex)[0];
    Real h21 = -4.0 * 100.0 * (*ex)[0];

    (*ehv)[0] = h11 * (*ev)[0] + h12 * (*ev)[1];
    (*ehv)[1] = h21 * (*ev)[0] + h22 * (*ev)[1];
  } 
#endif
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<std::vector<Real> > ehv
      = dynamic_cast<StdVector<Real>&>(hv).getVector();
    Ptr<const std::vector<Real> > ev
      = dynamic_cast<const StdVector<Real>&>(v).getVector();
    Ptr<const std::vector<Real> > ex
      = dynamic_cast<const StdVector<Real>&>(x).getVector();
    
    Real h11 = -4.0 * 100.0 * (*ex)[1] + 12.0 * 100.0 * std::pow((*ex)[0],2.0) + 2.0; 
    Real h22 = 2.0 * 100.0;
    Real h12 = -4.0 * 100.0 * (*ex)[0];
    Real h21 = -4.0 * 100.0 * (*ex)[0];

    (*ehv)[0] = 1.0/(h11*h22 - h12*h21) * (h22 * (*ev)[0] - h12 * (*ev)[1]);
    (*ehv)[1] = 1.0/(h11*h22 - h12*h21) * (-h21 * (*ev)[0] + h11 * (*ev)[1]);
  }
};

template<class Real>
class getHS1 : public TestProblem<Real> {
public:
  getHS1(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return makePtr<Objective_HS1<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem size
    int n = 2;
    // Get Initial Guess
    Ptr<std::vector<Real> > x0p = makePtr<std::vector<Real>>(n,0.0);
    (*x0p)[0] = -2.0; (*x0p)[1] = 1.0;
    return makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem size
    int n = 2;
    // Get Solution
    Ptr<std::vector<Real> > xp  = makePtr<std::vector<Real>>(n,0.0);
    (*xp)[0] = 1.0; (*xp)[1] = 1.0;
    return makePtr<StdVector<Real>>(xp);
  }

  Ptr<BoundConstraint<Real>> getBoundConstraint(void) const {
    // Problem size
    int n = 2;
    // Build lower bound
    Ptr<std::vector<Real> > lp = makePtr<std::vector<Real>>(n,0.0); 
    (*lp)[0] = ROL_NINF<Real>(); (*lp)[1] = -1.5;
    Ptr<Vector<Real> > l = makePtr<StdVector<Real>>(lp);
    // Build upper bound
    Ptr<std::vector<Real> > up = makePtr<std::vector<Real>>(n,0.0); 
    (*up)[0] = ROL_INF<Real>(); (*up)[1] = ROL_INF<Real>();
    Ptr<Vector<Real> > u = makePtr<StdVector<Real>>(up);
    // Instantiate BoundConstraint
    return makePtr<Bounds<Real>>(l,u);
  }
};

} // End ZOO Namespace
} // End ROL Namespace

#endif
