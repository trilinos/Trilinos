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
#include <iostream>

namespace ROL {

template<typename Real>
class PolyhedralProjection {
protected:
  const Ptr<BoundConstraint<Real>> bnd_;
  const Ptr<Constraint<Real>>      con_;
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

  const Ptr<Constraint<Real>> getLinearConstraint(void) const;

  const Ptr<BoundConstraint<Real>> getBoundConstraint(void) const;

  const Ptr<Vector<Real>> getMultiplier(void) const;

  const Ptr<Vector<Real>> getResidual(void) const;

}; // class PolyhedralProjection

} // namespace ROL

#include "ROL_PolyhedralProjection_Def.hpp"

#endif
