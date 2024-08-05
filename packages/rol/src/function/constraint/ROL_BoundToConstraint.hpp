// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BOUND_TO_CONSTRAINT_H
#define ROL_BOUND_TO_CONSTRAINT_H

#include "ROL_LowerBoundToConstraint.hpp"
#include "ROL_UpperBoundToConstraint.hpp"

/**  @ingroup func_group
     \class ROL::BoundToConstraint 
     \brief Provides an implementation for bound constraints.
*/

namespace ROL {

template<typename Real>
class BoundToConstraint : public Constraint<Real> { 
private:
  Ptr<Constraint<Real>> lo_;
  Ptr<Constraint<Real>> up_;
  Ptr<Vector<Real>> tmp_;

public:
  BoundToConstraint(BoundConstraint<Real> &bnd);
  BoundToConstraint(const Vector<Real> &lo, const Vector<Real> &up);

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                     const Vector<Real> &x, Real &tol) override;
};

}

#include "ROL_BoundToConstraint_Def.hpp"

#endif
