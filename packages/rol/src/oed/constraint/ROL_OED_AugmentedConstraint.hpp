// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_AUGMENTED_CONSTRAINT_HPP
#define ROL_OED_AUGMENTED_CONSTRAINT_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class AugmentedConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>> con_;
  const Ptr<Vector<Real>> rvec_;

public:
  // Probability Constraint
  AugmentedConstraint(const Vector<Real> &p, bool useScale = true, Real scale = Real(-1));

  //Budget Constraint
  AugmentedConstraint(const Ptr<Vector<Real>> &cost, Real budget);

  const Ptr<Vector<Real>> buildRangeVector() const { return rvec_->clone(); }

  void value(Vector<Real> &c,const Vector<Real> &x,Real &tol) override;
  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,
                     const Vector<Real> &x,Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,
                            const Vector<Real> &x,Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv,const Vector<Real> &u,
                           const Vector<Real> &v,const Vector<Real> &x,Real &tol) override;
};

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_AugmentedConstraint_Def.hpp"

#endif
