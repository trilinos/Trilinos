// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_H
#define ROL_AFFINE_HYPERPLANE_EQUALITY_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_Constraint.hpp"

#include <vector>
/** @ingroup func_group
    \class ROL::ScalarLinearConstraint
    \brief This equality constraint defines an affine hyperplane.

    ROL's scalar linear equality constraint interface implements
    \f[
       c(x) := \langle a, x\rangle_{\mathcal{X}^*,\mathcal{X}} - b = 0
    \f]
    where \f$a\in\mathcal{X}^*\f$ and \f$b\in\mathbb{R}\f$.  The range space of
    \f$c\f$ is an ROL::SingletonVector with dimension 1.

    Note: If \f$a\neq 0\f$ then there exists an explicit solution of the
    augmented system.  Namely,
    \f[
       v_1 = I^{-1}(b_1-av_2)
         \quad\text{and}\quad
       v_2 = \frac{(\langle a,I^{-1}b_1\rangle_{\mathcal{X}^*,\mathcal{X}}
               - b_2)}{\|a\|_{\mathcal{X}^*}^2}\,.
    \f]
    Moreover, note that \f$I^{-1}v\f$ for any \f$v\in\mathcal{X}^*\f$ is
    implemented in ROL as v.dual(). 

    ---
*/

namespace ROL {

template<typename Real>
class ScalarLinearConstraint : public Constraint<Real> {
private:
  const Ptr<const Vector<Real>> a_; ///< Dual vector defining hyperplane
  const Real b_;                    ///< Affine shift

public:
  ScalarLinearConstraint(const Ptr<const Vector<Real>> &a,
                         const Real b);

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v,
               const Vector<Real> &x,  Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                      const Vector<Real> &x,   Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u,
                     const Vector<Real> &v,    const Vector<Real> &x,
                           Real &tol) override;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                   const Vector<Real> &b1, const Vector<Real> &b2,
                                   const Vector<Real> &x,  Real &tol) override;

}; // class ScalarLinearConstraint

} // namespace ROL

#include "ROL_ScalarLinearConstraint_Def.hpp"

#endif
