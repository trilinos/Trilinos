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

#ifndef ROL_OBJECTIVE_H
#define ROL_OBJECTIVE_H

#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** \class ROL::Objective
    \brief Provides the interface to evaluate objective functions.
*/


namespace ROL {

template <class Real>
class Objective {
public:

  virtual ~Objective() {}

  /** \brief Update objective function.  
                x is a control, 
                flag = true if control is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  /** \brief Compute value.
  */
  virtual Real value( const Vector<Real> &x, Real &tol ) = 0;

  /** \brief Compute gradient.
  */
  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) ;

  /** \brief Compute directional derivative.
  */
  virtual Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) ;

  /** \brief Apply Hessian approximation to vector.
  */
  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol );

  /** \brief Apply inverse Hessian approximation to vector.
  */
  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {}

  /** \brief Apply preconditioner to vector.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v);
  }

  /** \brief Finite-difference gradient check.
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x,
                                                         const Vector<Real> &d,
                                                         const bool printToScreen = true,
                                                         const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

  /** \brief Finite-difference Hessian-applied-to-vector check.
  */
  virtual std::vector<std::vector<Real> > checkHessVec( const Vector<Real> &x,
                                                        const Vector<Real> &v,
                                                        const bool printToScreen = true,
                                                        const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

  /** \brief Hessian symmetry check.
  */
  virtual std::vector<Real> checkHessSym( const Vector<Real> &x,
                                          const Vector<Real> &v,
                                          const Vector<Real> &w,
                                          const bool printToScreen = true ) ;
  
  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

// include templated definitions
#include <ROL_ObjectiveDef.hpp>

#endif
