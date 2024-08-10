// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_I_HOM_OBJECTIVE_DEF_HPP
#define ROL_OED_I_HOM_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
I_Objective<Real>::I_Objective( const Ptr<BilinearConstraint<Real>> &con,
                                const Ptr<LinearObjective<Real>>    &obj,
                                const Ptr<Vector<Real>>             &state,
                                bool storage) {
  setConstraint(con);
  setObjective(obj);
  setStorage(storage);
  initialize(state);
}

template<typename Real>
Real I_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Get objective function value
  return getObjective()->value(state,z,tol);
}

template<typename Real>
void I_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Build gradient
  getConstraint()->applyAdjointJacobian_2(g,state,state,z,tol);
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void I_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Solve state sensitivity equation
  solve_state_sensitivity(v,z,tol);
  Vector<Real> &stateSens = *getStateSens();
  // Build hessVec
  getConstraint()->applyAdjointJacobian_2(hv,stateSens,state,z,tol);
  hv.scale(static_cast<Real>(-2));
}

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
