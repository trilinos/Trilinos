// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Rosenbrock's function.
    \author Created by D. Ridzal and D. Kouri.
 */

// Whether or not to use the exact Hessian-times-a-vector
#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_ROSENBROCK_HPP
#define ROL_ROSENBROCK_HPP

#include "ROL_StdVector.hpp"
#include "ROL_TestProblem.hpp"

namespace ROL {
namespace ZOO {

/** \brief Rosenbrock's function.
 */
template< class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real> >
class Objective_Rosenbrock : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;    

  typedef typename vector::size_type uint;

private:
  Real alpha_;

  Real const1_;
  Real const2_;

  template<class VectorType> 
  ROL::Ptr<const vector> getVector( const V& x ) {
    return dynamic_cast<const VectorType&>((x)).getVector();
  }

  template<class VectorType> 
  ROL::Ptr<vector> getVector( V& x ) {
    return dynamic_cast<VectorType&>(x).getVector();
  }

public:
  Objective_Rosenbrock(Real alpha = 100.0) : alpha_(alpha), const1_(100.0), const2_(20.0) {}

  Real value( const Vector<Real> &x, Real &tol ) {

    
    ROL::Ptr<const vector> xp = getVector<XPrim>(x);

    uint n = xp->size();
    Real val = 0;
    for( uint i=0; i<n/2; i++ ) {
      val += alpha_ * pow(pow((*xp)[2*i],2) - (*xp)[2*i+1], 2);
      val += pow((*xp)[2*i] - 1.0, 2);
    }

    //////  ADD INEXACTNESS
    //Real error = tol*(2.0*((Real)rand())/((Real)RAND_MAX)-1.0);
    //val += this->const1_*error; 

    return val;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    
    ROL::Ptr<const vector> xp = getVector<XPrim>(x);
    ROL::Ptr<vector> gp = getVector<XDual>(g);

    uint n = xp->size();
    for( uint i=0; i<n/2; i++ ) {
      (*gp)[2*i]   =  4.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1])*(*xp)[2*i] + 2.0*((*xp)[2*i]-1.0);
      (*gp)[2*i+1] = -2.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1]);

      //////  ADD INEXACTNESS
      //Real error0        = tol*(2.0*((Real)rand())/((Real)RAND_MAX)-1.0);
      //Real error1        = tol*(2.0*((Real)rand())/((Real)RAND_MAX)-1.0);
      //(*gp)[2*i]   += this->const2_*error0/std::sqrt(n);
      //(*gp)[2*i+1] += this->const2_*error1/std::sqrt(n);
    }
  }
#if USE_HESSVEC
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    
    ROL::Ptr<const vector> xp = getVector<XPrim>(x);
    ROL::Ptr<const vector> vp = getVector<XPrim>(v);
    ROL::Ptr<vector> hvp = getVector<XDual>(hv);

    uint n = xp->size();
    for( uint i=0; i<n/2; i++ ) {
      Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
      Real h12 = -4.0*alpha_*(*xp)[2*i];
      Real h22 = 2.0*alpha_;

      (*hvp)[2*i]   = h11*(*vp)[2*i] + h12*(*vp)[2*i+1];
      (*hvp)[2*i+1] = h12*(*vp)[2*i] + h22*(*vp)[2*i+1];
    }
  }
#endif
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  
     
  
    ROL::Ptr<const vector> xp = getVector<XPrim>(x);
    ROL::Ptr<const vector> vp = getVector<XDual>(v);
    ROL::Ptr<vector> hvp = getVector<XPrim>(hv);

    uint n = xp->size();
    for( uint i=0; i<n/2; i++ ) {
      Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
      Real h12 = -4.0*alpha_*(*xp)[2*i];
      Real h22 = 2.0*alpha_;

      (*hvp)[2*i]   = (1.0/(h11*h22-h12*h12))*( h22*(*vp)[2*i] - h12*(*vp)[2*i+1]);
      (*hvp)[2*i+1] = (1.0/(h11*h22-h12*h12))*(-h12*(*vp)[2*i] + h11*(*vp)[2*i+1]);
    }
  }
};

template<class Real>
class getRosenbrock : public TestProblem<Real> {
public:
  getRosenbrock(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Instantiate Objective Function
    return ROL::makePtr<Objective_Rosenbrock<Real>>();
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 100;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n/2; i++ ) {
      (*x0p)[2*i]   = -1.2;
      (*x0p)[2*i+1] =  1.0;
    }
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 100;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) {
      (*xp)[i] = 1.0;
    }
    return ROL::makePtr<StdVector<Real>>(xp);
  }
};

}// End ZOO Namespace
}// End ROL Namespace

#endif
