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

#ifndef ROL_INEQUALITY_CONSTRAINT_H
#define ROL_INEQUALITY_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** \class ROL::InequalityConstraint
    \brief Provides the interface to evaluate an inequality constraint function.
*/


namespace ROL {

template <class Real>
class InequalityConstraint {
private:
  bool activated_;

public:

  virtual ~InequalityConstraint() {}

  InequalityConstraint(void) : activated_(true) {}

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if control is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  /** \brief Project optimization variables onto constraint set.
                x is the optimization variable
  */
  virtual void project( Vector<Real> &x ) {}

  /** \brief Remove active set variables that are also in the binding set.
                v is the vector to be pruned 
                g is the gradient of the objective function at x
                x is the optimization variable
                eps is the active set tolerance
  */
  virtual void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {}

  /** \brief Remove active set variables.
                v is the vector to be pruned 
                x is the optimization variable
                eps is the active set tolerance
  */
  virtual void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {}

  /** \brief Check if the vector, v, is feasible
  */
  virtual bool isFeasible( const Vector<Real> &v ) { return true; }

  /** \brief Evaluate constraint.
  */
  virtual void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {}
 
  /** \brief Apply constraint Jacobian or its adjoint.
  */
  virtual void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, const bool &adj, Real &tol );

  /** \brief Computes the action of the operator W that is onto
             the null space (kernel) of the contraint Jacobian.
  */
  virtual void maptoJacobianKernel( Vector<Real> &wv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {}

  /** \brief Apply constraint Hessian to (v,u): c''(x)(v,u) = (c''(x)u)v.
  */
  virtual void applyHessian( Vector<Real> &huv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol );

  /** \brief Turn on constraints 
  */
  void activate(void)    { this->activated_ = true;  }

  /** \brief Turn off constraints
  */
  void deactivate(void)  { this->activated_ = false; }

  /** \brief Check if constraints are on
  */
  bool isActivated(void) { return this->activated_;  }

}; // class InequalityConstraint

} // namespace ROL

#include "ROL_InequalityConstraintDef.hpp"

#endif
