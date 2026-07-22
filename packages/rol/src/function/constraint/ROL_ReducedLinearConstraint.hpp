// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCED_LINEAR_CONSTRAINT_H
#define ROL_REDUCED_LINEAR_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::ReducedLinearConstraint
    \brief Reduce the input of a linear constraint based on the active set
	   associated with a vector \f$x\f$, i.e., let \f$\mathcal{I}\f$ denote
	   the inactive set associated with \f$x\f$ and the bounds
           \f$\ell\le u\f$, then

           \f[ C(v) = c(v_\mathcal{I}), \f]

           where \f$v_\mathcal{I}\f$ denotes the vector that is equal to
           \f$v\f$ on \f$\mathcal{I}\f$ and zero otherwise.

*/

namespace ROL {

template<typename Real>
class ReducedLinearConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>>      con_;
  const Ptr<BoundConstraint<Real>> bnd_;
  Ptr<const Vector<Real>>          x_;
  const Ptr<Vector<Real>>          prim_;

public:
  ReducedLinearConstraint(const Ptr<Constraint<Real>> &con,
                          const Ptr<BoundConstraint<Real>> &bnd,
                          const Ptr<const Vector<Real>> &x);

  void setX(const Ptr<const Vector<Real>> &x);

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) override;
}; // class ReducedLinearConstraint

} // namespace ROL

#include "ROL_ReducedLinearConstraint_Def.hpp"

#endif // ROL_REDUCED_LINEAR_CONSTRAINT_H
