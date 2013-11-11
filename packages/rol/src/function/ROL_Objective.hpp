//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_OBJECTIVE_H
#define ROL_OBJECTIVE_H

#include "ROL_Vector.hpp"
#include <iostream>

/** \class ROL::Objective
    \brief Provides the interface to evaluate objective functions.
*/


namespace ROL {

template <class Real>
class Objective {
public:

  virtual ~Objective() {}

  /** \brief Compute value.
  */
  virtual Real value( const Vector<Real> &x ) = 0;

  /** \brief Compute gradient.
  */
  virtual void gradient( Vector<Real> &g, const Vector<Real> &x ) ;

  /** \brief Compute directional derivative.
  */
  virtual Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, const Real eta = 1.0 ) ;

  /** \brief Apply Hessian approximation to vector.
  */
  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x );

  /** \brief Apply inverse Hessian approximation to vector.
  */
  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x ) {}

  /** \brief Apply preconditioner to vector.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x ) {
    Pv.set(v);
  }

  /** \brief Finite-difference gradient check.
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x, const Vector<Real> &d, const bool printToScreen = true ) ;
  
  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

// include templated definitions
#include <ROL_ObjectiveDef.hpp>

#endif
