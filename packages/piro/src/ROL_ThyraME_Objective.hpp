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

#ifndef ROL_THYRAME_OBJECTIVE_H
#define ROL_THYRAME_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::Objective
    \brief Provides the interface to evaluate objective functions.

    ROL's objective function interface is designed for Fr$eacute;chet differentiable 
    functionals \f$f:\mathcal{X}\to\mathbb{R}\f$, where \f$\mathcal{X}\f$ is a Banach
    space.  The basic operator interace, to be implemented by the user, requires:
    \li #value -- objective function evaluation.

    It is strongly recommended that the user additionally overload:
    \li #gradient -- the objective function gradient -- the default is a finite-difference approximation;
    \li #hessVec  -- the action of the Hessian -- the default is a finite-difference approximation.

    The user may also overload:
    \li #update     -- update the objective function at each new iteration;
    \li #dirDeriv   -- compute the directional derivative -- the default is a finite-difference approximation;
    \li #invHessVec -- the action of the inverse Hessian;
    \li #precond    -- the action of a preconditioner for the Hessian.

    ---
*/


namespace ROL {

template <class Real>
class ThyraME_Objective : public Objective<Real> {
public:

  ThyraME_Objective() {};

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    return 0.5 * x.dot(x);
  };

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
     g.set(x);
  };

}; // class Objective

} // namespace ROL
#endif
