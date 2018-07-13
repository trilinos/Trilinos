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
