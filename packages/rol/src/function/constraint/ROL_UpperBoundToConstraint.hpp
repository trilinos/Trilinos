// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UPPER_BOUND_TO_CONSTRAINT_H
#define ROL_UPPER_BOUND_TO_CONSTRAINT_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"

/**  @ingroup func_group
     \class ROL::UpperBoundToConstraint 
     \brief Provides an implementation for upper bound constraints.
*/

namespace ROL {

template<typename Real>
class UpperBoundToConstraint : public Constraint<Real> { 
private:
  Ptr<Vector<Real>> up_;

public:
  UpperBoundToConstraint(BoundConstraint<Real> &bnd);
  UpperBoundToConstraint(const Vector<Real> &up);

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) override;
};

}

#include "ROL_UpperBoundToConstraint_Def.hpp"

#endif
