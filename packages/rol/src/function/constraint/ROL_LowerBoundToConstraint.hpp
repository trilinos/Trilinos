// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LOWER_BOUND_TO_CONSTRAINT_H
#define ROL_LOWER_BOUND_TO_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"

/**  @ingroup func_group
     \class ROL::LowerBoundToConstraint 
     \brief Provides an implementation for lower bound constraints.
*/

namespace ROL {

template<typename Real>
class LowerBoundToConstraint : public Constraint<Real> {
private:
  Ptr<Vector<Real>> lo_;

public:
  LowerBoundToConstraint(BoundConstraint<Real> &bnd);
  LowerBoundToConstraint(const Vector<Real> &lo);

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) override;
};

}

#include "ROL_LowerBoundToConstraint_Def.hpp"

#endif
