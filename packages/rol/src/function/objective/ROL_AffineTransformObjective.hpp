// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_TRANSFORM_OBJECTIVE_H
#define ROL_AFFINE_TRANSFORM_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_LinearConstraint.hpp"
#include "ROL_VectorController.hpp"

/** @ingroup func_group
    \class ROL::AffineTransformObjective
    \brief Compose an objective function with an affine transformation, i.e.,

    \f[ F(x) = f(Ax+b). \f]

*/

namespace ROL {

template<typename Real>
class AffineTransformObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>>   obj_;
  const Ptr<Constraint<Real>> acon_;

  Ptr<VectorController<Real>> storage_;
  Ptr<Vector<Real>> primal_, dual_, Av_;

public:
  virtual ~AffineTransformObjective() {}
  AffineTransformObjective(const Ptr<Objective<Real>>        &obj,
                           const Ptr<Constraint<Real>>       &acon,
                           const Vector<Real>                &range,
                           const Ptr<VectorController<Real>> &storage = nullPtr);
  AffineTransformObjective(const Ptr<Objective<Real>>        &obj,
                           const Ptr<LinearConstraint<Real>> &acon,
                           const Ptr<VectorController<Real>> &storage = nullPtr);
  AffineTransformObjective(const Ptr<Objective<Real>>            &obj,
                           const Ptr<const LinearOperator<Real>> &A,
                           const Ptr<const Vector<Real>>         &b,
                           const Ptr<VectorController<Real>>     &storage = nullPtr);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

public:
  void setParameter(const std::vector<Real> &param) override;

private:
  Ptr<const Vector<Real>> transform(const Vector<Real> &x);

}; // class AffineTransformObjective

} // namespace ROL

#include "ROL_AffineTransformObjective_Def.hpp"

#endif // ROL_AFFINE_TRANSFORM_OBJECTIVE_H
