// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDEQUALITY_CONSTRAINT_H
#define ROL_STDEQUALITY_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::StdConstraint
    \brief Defines the equality constraint operator interface for StdVectors

*/

namespace ROL {

template<typename Real>
class StdConstraint : public virtual Constraint<Real> {
public:
  virtual ~StdConstraint() {}

  using Constraint<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}  
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  virtual void update( const std::vector<Real> &x, UpdateType type, int iter = -1 ) {}

  using Constraint<Real>::value;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  virtual void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) = 0;

  using Constraint<Real>::applyJacobian;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                             const Vector<Real> &x, Real &tol) override;
  virtual void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                              const std::vector<Real> &x, Real &tol );

  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                                    const Vector<Real> &x, Real &tol) override;
   virtual void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                                      const std::vector<Real> &x, Real &tol );

  using Constraint<Real>::applyAdjointHessian;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) override;
  virtual void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                                    const std::vector<Real> &v, const std::vector<Real> &x,
                                    Real &tol );

  using Constraint<Real>::solveAugmentedSystem;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                         const Vector<Real> &b1, const Vector<Real> &b2,
                                         const Vector<Real> &x, Real &tol) override;
  virtual std::vector<Real> solveAugmentedSystem( std::vector<Real> &v1, std::vector<Real> &v2,
                                                  const std::vector<Real> &b1, const std::vector<Real> &b2,
                                                  const std::vector<Real> &x, Real tol );

  using Constraint<Real>::applyPreconditioner;
  void applyPreconditioner(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x,
                           const Vector<Real> &g, Real &tol) override;
  virtual void applyPreconditioner( std::vector<Real> &pv, const std::vector<Real> &v,
                                    const std::vector<Real> &x, const std::vector<Real> &g, Real &tol );

}; // class StdConstraint

} // namespace ROL

#include "ROL_StdConstraint_Def.hpp"

#endif
