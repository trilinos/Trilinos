// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_TRANSFORM_CONSTRAINT_H
#define ROL_AFFINE_TRANSFORM_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_LinearConstraint.hpp"
#include "ROL_VectorController.hpp"

/** @ingroup func_group
    \class ROL::AffineTransformConstraint
    \brief Compose a constraint operator with an affine transformation, i.e.,

    \f[ C(x) = c(Ax+b). \f]

*/

namespace ROL {

template <class Real>
class AffineTransformConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>>  con_;
  const Ptr<Constraint<Real>> acon_;

  Ptr<VectorController<Real>> storage_;
  Ptr<Vector<Real>> primal_, dual_, Av_;

public:
  virtual ~AffineTransformConstraint() {}
  AffineTransformConstraint(const Ptr<Constraint<Real>>       &con,
                            const Ptr<Constraint<Real>>       &acon,
                            const Vector<Real>                &range,
                            const Ptr<VectorController<Real>> &storage = nullPtr);
  AffineTransformConstraint(const Ptr<Constraint<Real>>       &con,
                            const Ptr<LinearConstraint<Real>> &acon,
                            const Ptr<VectorController<Real>> &storage = nullPtr);
  AffineTransformConstraint(const Ptr<Constraint<Real>>           &con,
                            const Ptr<const LinearOperator<Real>> &A,
                            const Ptr<const Vector<Real>>         &b,
                            const Ptr<VectorController<Real>>     &storage = nullPtr);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) override;
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

private:
  Ptr<const Vector<Real>> transform(const Vector<Real> &x);

}; // class AffineTransformConstraint

} // namespace ROL

#include "ROL_AffineTransformConstraint_Def.hpp"

#endif // ROL_AFFINE_TRANSFORM_OBJECTIVE_H
