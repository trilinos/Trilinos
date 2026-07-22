// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKLESS_CONSTRAINT_H
#define ROL_RISKLESS_CONSTRAINT_H

#include "ROL_RiskVector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

template <class Real>
class RiskLessConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>> con_;

public:
  RiskLessConstraint(const Ptr<Constraint<Real>> &con);

  void update(const Vector<Real> &x, UpdateType type, int iter = -1) override;
  void update(const Vector<Real> &x, bool flag = true, int iter = -1) override;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;

public:
  void setParameter(const std::vector<Real> &param) override;

}; // class RiskLessConstraint

} // namespace ROL

#include "ROL_RiskLessConstraint_Def.hpp"

#endif
