// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ELASTICLINEARCONSTRAINT_H
#define ROL_ELASTICLINEARCONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::ElasticLinearConstraint
    \brief Defines the general affine constraint with the form \f$c(x)=g(x) + g'(x)s + u - v\f$.

    ---
*/

namespace ROL {

template<typename Real>
class ElasticLinearConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>> con_;
  const Ptr<Vector<Real>> x_, c_, tmp_;

public:
  ElasticLinearConstraint(const Ptr<const Vector<Real>> &x,
                          const Ptr<Constraint<Real>>   &con,
                          const Ptr<const Vector<Real>> &c);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;

  void setAnchor(const Ptr<const Vector<Real>> &x);

}; // class ElasticLinearConstraint

} // namespace ROL

#include "ROL_ElasticLinearConstraint_Def.hpp"

#endif
