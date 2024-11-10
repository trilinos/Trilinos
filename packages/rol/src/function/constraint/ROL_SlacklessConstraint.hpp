// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SLACKLESSCONSTRAINT_HPP
#define ROL_SLACKLESSCONSTRAINT_HPP

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
 *  \class ROL::SlacklessConstraint
 *  \brief This class strips out the slack variables from constraint evaluations
 *         to create the new constraint  \f$ C(x,s) = c(x) \f$
 */

namespace ROL {

template<typename Real> 
class SlacklessConstraint : public Constraint<Real> {
private: 
  const Ptr<Constraint<Real>> con_;

public:
  SlacklessConstraint( const Ptr<Constraint<Real>> &con );
 
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) override;
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &dualv, Real &tol ) override;
  void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

// Definitions for parametrized (stochastic) constraint functions
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  Ptr<Vector<Real>> getOpt( Vector<Real> &xs ) const;
  Ptr<const Vector<Real>> getOpt( const Vector<Real> &xs ) const;
  void zeroSlack( Vector<Real> &x ) const;

}; // class SlacklessConstraint 

} // namespace ROL

#include "ROL_SlacklessConstraint_Def.hpp"

#endif // ROL__SLACKLESSCONSTRAINT_HPP

