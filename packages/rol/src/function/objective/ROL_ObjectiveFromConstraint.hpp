// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_FROM_CONSTRAINT_H
#define ROL_OBJECTIVE_FROM_CONSTRAINT_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ObjectiveFromConstraint
    \brief Form an objective function from a ROL::Constraint and a
           vector in the dual constraint space \f$\lambda\in \mathcal{C}^\ast\f$

    \f[ f(x;\lambda) = \langle \lambda, c(x)\rangle_{\mathcal{C}^*,\mathcal{C}} \f]
*/

namespace ROL {

template<typename Real>
class ObjectiveFromConstraint : public Objective<Real> {
private:
  Ptr<Constraint<Real>> con_;
  Ptr<Vector<Real>>     l_;      // Lagrange multiplier 
  Ptr<Vector<Real>>     c_;      // Constraint vector

public:
  virtual ~ObjectiveFromConstraint() {}
  ObjectiveFromConstraint( const Ptr<Constraint<Real>> &con, 
                           const Vector<Real> &l );

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void updateMultiplier( const Vector<Real> &l );

}; // class ObjectiveFromConstraint

} // namespace ROL

#include "ROL_ObjectiveFromConstraint_Def.hpp"

#endif // ROL_OBJECTIVE_FROM_CONSTRAINT_H
