// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PQNOBJECTIVE_H
#define ROL_PQNOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Secant.hpp"

/** @ingroup func_group
    \class ROL::PQNObjective
    \brief Provides the interface to evaluate the quadratic quasi-Newton objective.

    This class implements the PQN quasi-Newton objective for use with
    ROL::TypeB::QuasiNewtonAlgorithm.  Given a function
    \f$f:\mathcal{X}\to\mathbb{R}\f$ and a Hessian approximation
    \f$B_k:\mathcal{X}}\to\mathcal{X}^*\f$, the functional is
    \f[
       q_k(x) = \frac{1}{2}\langle B_k (x-x_k),(x-x_k)\rangle_{\mathcal{X}^*,\mathcal{X}}
                + \langle f'(x_k), (x-x_k)\rangle_{\mathcal{X}^*,\mathcal{X}}.
    \f]

    ---
*/


namespace ROL {

template<typename Real>
class PQNObjective : public Objective<Real> {
private:
  const Ptr<Secant<Real>> secant_;
  const Ptr<Vector<Real>> x_, g_, pwa_, dwa_;

public:
  PQNObjective(const Ptr<Secant<Real>> &secant,
               const Vector<Real> &x,
               const Vector<Real> &g);

  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

  void setAnchor(const Vector<Real> &x, const Vector<Real> &g);
}; // class PQNObjective

} // namespace ROL

#include "ROL_PQNObjective_Def.hpp"

#endif
