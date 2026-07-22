// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATIC_OBJECTIVE_H
#define ROL_QUADRATIC_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::QuadraticObjective
    \brief Provides the interface to evaluate quadratic objective functions.

    This class implements the quadratic objective function
    \f[
       f(x) = \frac{1}{2}\langle Hx, x\rangle_{\mathcal{X}^*,\mathcal{X}}
            + \langle g,  x\rangle_{\mathcal{X}^*,\mathcal{X}}
            + c
    \f]
    for fixed \f$H\in\mathcal{L}(\mathcal{X},\mathcal{X}^*)\f$,
    \f$g\in\mathcal{X}^*\f$, and \f$c\in\mathbb{R}\f$.

    ---
*/

namespace ROL {

template<typename Real>
class QuadraticObjective : public Objective<Real> {
private:
  const Ptr<const LinearOperator<Real>> H_;
  const Ptr<const Vector<Real>> g_;
  const Real c_;
  Ptr<Vector<Real>> tmp_;

public:
  QuadraticObjective(const Ptr<const LinearOperator<Real>> &H,
                     const Ptr<const Vector<Real>>         &g,
                     Real                                   c = Real(0));

  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

}; // class QuadraticObjective

} // namespace ROL

#include "ROL_QuadraticObjective_Def.hpp"

#endif
