// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINTFROMOBJECTIVE_H
#define ROL_CONSTRAINTFROMOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_SingletonVector.hpp"

/** @ingroup func_group
    \class ROL::ConstraintFromObjective
    \brief Creates a constraint from an objective function and a offset value

    Example:  Suppose we have an objective function f(x) and we wish to impose,
	      e.g., a condition f(x)-offset = 0, then this class creates the
              scalar constraint c(x) = f(x)-offset 
*/


namespace ROL {

template<typename Real> 
class ConstraintFromObjective : public Constraint<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>>          dualVector_;
  const Real                 offset_;
  bool                       isDualInitialized_;

public:
  ConstraintFromObjective( const Ptr<Objective<Real>> &obj, const Real offset = 0 );

  const Ptr<Objective<Real>> getObjective(void) const;

  void setParameter( const std::vector<Real> &param ) override;

  void update( const Vector<Real>& x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real>& x, bool flag = true, int iter = -1 ) override;
  void value( Vector<Real>& c, const Vector<Real>& x, Real& tol ) override;
  void applyJacobian( Vector<Real>& jv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;
  void applyAdjointJacobian( Vector<Real>& ajv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;
  void applyAdjointHessian( Vector<Real>& ahuv, const Vector<Real>& u, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;

private:
  Real getValue( const Vector<Real>& x ); 
  void setValue( Vector<Real>& x, Real val );

}; // ConstraintFromObjective

} // namespace ROL

#include "ROL_ConstraintFromObjective_Def.hpp"

#endif // ROL_CONSTRAINTFROMOBJECTIVE_H
