// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEARCONSTRAINT_H
#define ROL_LINEARCONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::LinearConstraint
    \brief Defines the general affine constraint with the form \f$c(x)=Ax+b\f$.

    ---
*/

namespace ROL {

template<typename Real>
class LinearConstraint : public Constraint<Real> {
private:
  const Ptr<const LinearOperator<Real>> A_;
  const Ptr<const Vector<Real>> b_;

public:
  LinearConstraint(const Ptr<const LinearOperator<Real>> &A,
                   const Ptr<const Vector<Real>> &b);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;

  Ptr<Vector<Real>> createRangeSpaceVector(void) const;

}; // class LinearConstraint

} // namespace ROL

#include "ROL_LinearConstraint_Def.hpp"

#endif
