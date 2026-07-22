// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAROPERATOR_H
#define ROL_LINEAROPERATOR_H

#include "ROL_Vector.hpp"

/** @ingroup func_group
    \class ROL::LinearOperator
    \brief Provides the interface to apply a linear operator.

    ROL's linear operator interface is designed to interface with ROL's Krylov methods.
    These linear operators often represent projected Hessians or preconditioners.  
    The basic operator interace, to be implemented by the user, requires:
    \li #apply -- apply operator to a vector.

    The user may also implement:
    \li #update -- update the state of the linear operator.
    \li #applyInverse -- apply the inverse operator to a vector.
    \li #applyAdjoint -- apply the adjoint of the operator to a vector.
    \li #applyAdjointInverse -- apply the adjoint of the inverse operator to a vector.

    ---
*/


namespace ROL {

template <class Real>
class LinearOperator {
public:

  virtual ~LinearOperator() {}

  /** \brief Update linear operator. 

      This function updates the linear operator at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  /** \brief Apply linear operator.

      This function applies the linear operator to a vector.
      @param[out]         Hv  is the output vector.
      @param[in]          v   is the input vector.
      @param[in]          tol is a tolerance for inexact linear operator application.
  */
  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const = 0;

  /** \brief Apply inverse of linear operator.

      This function applies the inverse of linear operator to a vector.
      @param[out]         Hv  is the output vector.
      @param[in]          v   is the input vector.
      @param[in]          tol is a tolerance for inexact linear operator application.
  */
  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v);
  }

  /** \brief Apply adjoint of linear operator.

      This function applies the adjoint of a linear operator to a vector. Default behavior
      assumes operator is self-adjoint.

      @param[out]         Hv  is the output vector.
      @param[in]          v   is the input vector.
      @param[in]          tol is a tolerance for inexact linear operator application.
  */
  virtual void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    apply(Hv,v,tol);
  }

  /** \brief Apply adjoint of the inverse linear operator.

      This function applies the adjoint of the inverse linear operator to a vector.
      Default behavior assumes operator is self-adjoint.
      @param[out]         Hv  is the output vector.
      @param[in]          v   is the input vector.
      @param[in]          tol is a tolerance for inexact linear operator application. 

  */
  virtual void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    applyInverse(Hv,v,tol);
  }

}; // class LinearOperator

} // namespace ROL

#endif
