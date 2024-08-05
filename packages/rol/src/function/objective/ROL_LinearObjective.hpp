// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAR_OBJECTIVE_H
#define ROL_LINEAR_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Ptr.hpp"

/** @ingroup func_group
    \class ROL::LinearObjective
    \brief Provides the interface to evaluate linear objective functions.

    This class implements the linear objective function
    \f[
       f(x) = \langle c, x\rangle_{\mathcal{X}^*,\mathcal{X}}
    \f]
    for fixed \f$c\in\mathcal{X}^*\f$.

    ---
*/

namespace ROL {

template<typename Real>
class LinearObjective : public Objective<Real> {
private:
  const Ptr<const Vector<Real>> cost_, dual_cost_;

public:
  LinearObjective(const Ptr<const Vector<Real>> &cost);

  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

}; // class LinearObjective

} // namespace ROL

#include "ROL_LinearObjective_Def.hpp"

#endif
