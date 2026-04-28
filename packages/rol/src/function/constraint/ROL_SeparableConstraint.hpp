// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SEPARABLE_CONSTRAINT_H
#define ROL_SEPARABLE_CONSTRAINT_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Ptr.hpp"
#include <vector>

namespace ROL {

/** @ingroup func_group
 *  \class ROL::SeparableConstraint
 *  \brief Array of constraints with entries that depend on indepdent subvectors.
 *
 */

template<typename Real>
class SeparableConstraint : public Constraint<Real> {
private:
  const std::vector<Ptr<Constraint<Real>>> con_;
  const unsigned size_;

public:
  SeparableConstraint(const std::vector<Ptr<Constraint<Real>>> &cvec);
  Ptr<Constraint<Real>> get(unsigned ind = 0) const;
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) override;
  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) override;
  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) override;
  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) override;
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) override;
  void setParameter(const std::vector<Real> &param) override;

}; // class SeparableConstraint

} // namespace ROL

#include "ROL_SeparableConstraint_Def.hpp"

#endif
