// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_C_HOM_OBJECTIVE_DEF_HPP
#define ROL_OED_C_HOM_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
C_Objective<Real>::C_Objective( const Ptr<BilinearConstraint<Real>> &con,
             const Ptr<LinearObjective<Real>>  &obj,
             const Ptr<Vector<Real>>           &state,
             const bool storage) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real C_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(pnull_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Get objective function value
  return ObjectiveBase<Real,std::vector<Real>>::getObjective()->value(state,z,tol);
}

template<typename Real>
void C_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(pnull_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Build gradient
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(g,state,state,z,tol);
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void C_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(pnull_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Solve state sensitivity equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity(v,z,tol);
  Vector<Real> &stateSens = *ObjectiveBase<Real,std::vector<Real>>::getStateSens();
  // Build hessVec
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(hv,stateSens,state,z,tol);
  hv.scale(static_cast<Real>(-2));
}

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
