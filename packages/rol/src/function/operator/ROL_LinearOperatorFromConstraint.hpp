// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
