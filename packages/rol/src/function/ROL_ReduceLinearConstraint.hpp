// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCE_LINEAR_CONSTRAINT_H
#define ROL_REDUCE_LINEAR_CONSTRAINT_H

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_AffineTransformConstraint.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_RangeSpaceOperator.hpp"

/** @ingroup func_group
    \class ROL::ReduceLinearConstraint
    \brief Performs null-space transformation for reducible linear equality
           constraints.

    ---
*/

namespace ROL {

template<typename Real>
class ReduceLinearConstraint {
private:
  const Ptr<Constraint<Real>> lcon_;
  const Ptr<Vector<Real>>     x_;
  const Ptr<VectorController<Real>>  storage_;
  const Ptr<NullSpaceOperator<Real>> nsop_;

public:
  virtual ~ReduceLinearConstraint(void) {}

  ReduceLinearConstraint(const Ptr<Constraint<Real>>   &lcon,
                         const Ptr<Vector<Real>>       &x,
                         const Ptr<const Vector<Real>> &c);

  Ptr<Objective<Real>> transform(const Ptr<Objective<Real>> &obj) const;
  Ptr<Constraint<Real>> transform(const Ptr<Constraint<Real>> &con) const;
  Ptr<Constraint<Real>> getLinearConstraint(void) const;
  Ptr<const Vector<Real>> getFeasibleVector(void) const;
  void project(Vector<Real> &x, const Vector<Real> &y) const;
  void project(const Ptr<Vector<Real>> &x, const Ptr<const Vector<Real>> &y) const;

private:
  void feasible(const Ptr<const Vector<Real>> &c);

}; // class ReduceLinearConstraint

} // namespace ROL

#include "ROL_ReduceLinearConstraint_Def.hpp"

#endif
