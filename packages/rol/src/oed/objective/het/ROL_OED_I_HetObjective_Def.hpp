// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_I_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_I_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
I_Objective<Real>::I_Objective( const Ptr<BilinearConstraint<Real>>   &con,
             const Ptr<QuadraticObjective<Real>> &obj,
             const Ptr<Vector<Real>>        &state,
             bool storage) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real I_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  // Solve state equation
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Get objective function value
  return getObjective()->value(state,z,tol);
}

template<typename Real>
void I_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (dp_==nullPtr) dp_ = g.clone();
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  // Solve state equation
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Solve adjoint equation
  solve_adjoint_equation(param,z,tol);
  Vector<Real> &adjoint = *getAdjoint();
  // Build gradient
  getObjective()->gradient_2(*dp_,state,z,tol);
  getConstraint()->applyAdjointJacobian_2(g,adjoint,state,z,tol);
  g.plus(*dp_);
}

template<typename Real>
void I_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (dp_==nullPtr) dp_ = hv.clone();
  std::vector<Real> param = ObjectiveBase<Real,std::vector<Real>>::getParameter();
  // Solve state equation
  solve_state_equation(param,z,tol);
  Vector<Real> &state = *getState();
  // Solve adjoint equation
  solve_adjoint_equation(param,z,tol);
  Vector<Real> &adjoint = *getAdjoint();
  // Solve state sensitivity equation
  solve_state_sensitivity(v,z,tol);
  Vector<Real> &stateSens = *getStateSens();
  // Solve adjoint sensitivity equation
  solve_adjoint_sensitivity(v,z,tol);
  Vector<Real> &adjointSens = *getAdjointSens();
  // Build hessVec
  getConstraint()->applyAdjointJacobian_2(hv,adjointSens,state,z,tol);
  getObjective()->hessVec_21(*dp_,stateSens,state,z,tol);
  hv.plus(*dp_);
  getConstraint()->applyAdjointHessian_12(*dp_,adjoint,stateSens,state,z,tol);
  hv.plus(*dp_);
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
