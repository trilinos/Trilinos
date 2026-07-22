// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_STATE_H
#define ROL_CONSTRAINT_STATE_H

#include "ROL_Constraint_SimOpt.hpp"

namespace ROL {

template<typename Real>
class SimConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint_SimOpt<Real>> con_;
  const Ptr<const Vector<Real>> z_;
  const bool inSolve_;
  Ptr<Vector<Real>> ijv_;
  bool init_;

public:
  SimConstraint(const Ptr<Constraint_SimOpt<Real>> &con,
                const Ptr<const Vector<Real>> &z,
                bool inSolve = false);

  void update( const Vector<Real> &u, bool flag = true, int iter = -1 ) override;
  void update( const Vector<Real> &u, UpdateType type, int iter = -1 ) override;
  void value(Vector<Real> &c,const Vector<Real> &u,Real &tol) override;
  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) override;
  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,const Vector<Real> &u,Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,const Vector<Real> &u,Real &tol) override;
  void applyPreconditioner(Vector<Real> &pv,const Vector<Real> &v,const Vector<Real> &u,const Vector<Real> &g,Real &tol) override;

  // Definitions for parametrized (stochastic) equality constraints
  void setParameter(const std::vector<Real> &param) override;

}; // class SimConstraint

} // namespace ROL

#include "ROL_SimConstraint_Def.hpp"

#endif
