// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_FSSOLVER_H
#define ROL_OBJECTIVE_FSSOLVER_H

#include "ROL_Objective.hpp"

namespace ROL {

template<typename Real>
class Objective_FSsolver : public Objective<Real> {
public:
  Real value( const Vector<Real> &u, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &u, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, Real &tol ) override;

}; // class Objective_FSsolver

} // namespace ROL

#include "ROL_Objective_FSsolver_Def.hpp"

#endif
