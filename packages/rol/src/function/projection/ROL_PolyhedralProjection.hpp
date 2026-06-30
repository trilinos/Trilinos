// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_POLYHEDRALPROJECTION_H
#define ROL_POLYHEDRALPROJECTION_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_Projection.hpp"
#include "ROL_ReducedLinearConstraint.hpp"
#include <iostream>

namespace ROL {

template<typename Real>
class PolyhedralProjection : public Projection<Real> {
protected:
  const Ptr<BoundConstraint<Real>>   bnd_;
  const Ptr<Constraint<Real>>        con_;
  Ptr<ReducedLinearConstraint<Real>> rcon_;  // con_ restricted to current active variables
  Ptr<NullSpaceOperator<Real>>       ns_;    // null space projection onto reduced equality constraint
  Ptr<Vector<Real>> xprim_, xdual_, mul_, res_;

public:
  virtual ~PolyhedralProjection() {}

  PolyhedralProjection(const Ptr<BoundConstraint<Real>> &bnd);

  PolyhedralProjection(const Vector<Real>               &xprim,
                       const Vector<Real>               &xdual,
                       const Ptr<BoundConstraint<Real>> &bnd,
                       const Ptr<Constraint<Real>>      &con,
                       const Vector<Real>               &mul,
                       const Vector<Real>               &res);

  virtual void project(Vector<Real> &x, std::ostream &stream = std::cout);

  virtual void applyJacobian(Vector<Real> &v, const Vector<Real> &x);

  const Ptr<Constraint<Real>> getLinearConstraint(void) const;

  const Ptr<BoundConstraint<Real>> getBoundConstraint(void) const;

  const Ptr<Vector<Real>> getMultiplier(void) const;

  const Ptr<Vector<Real>> getResidual(void) const;

}; // class PolyhedralProjection

} // namespace ROL

#include "ROL_PolyhedralProjection_Def.hpp"

#endif
