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

#ifndef ROL_LINEAROPERATOR_FROM_EQUALITYCONSTRAINT_H
#define ROL_LINEAROPERATOR_FROM_EQUALITYCONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::LinearOperatorFromConstraint
    \brief A simple wrapper which allows application of constraint Jacobians 
           through the LinearOperator interface

    ---
*/


namespace ROL {

template <class Real>
class LinearOperatorFromConstraint : public LinearOperator<Real> {
private:
  const ROL::Ptr<const Vector<Real> > x_;
  ROL::Ptr<Constraint<Real> > con_;
 

public:

  LinearOperatorFromConstraint( const ROL::Ptr<const Vector<Real> > &x, 
                                        const ROL::Ptr<Constraint<Real> > &con ) : 
                                        x_(x), con_(con) {
  }

  virtual ~LinearOperatorFromConstraint() {}

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    con_->applyJacobian(Hv,v,*x_,tol);
  }

  // Not implemented for generic equality constraint
  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v);
  }

}; // class LinearOperator

} // namespace ROL

#endif
