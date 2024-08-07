// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_C_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_C_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
C_Objective<Real>::C_Objective( const Ptr<BilinearConstraint<Real>>   &con,
             const Ptr<QuadraticObjective<Real>> &obj,
             const Ptr<Vector<Real>>             &state,
             bool storage) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real C_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Get objective function value
  return ObjectiveBase<Real,std::vector<Real>>::getObjective()->value(state,z,tol);
}

template<typename Real>
void C_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (dp_ == nullPtr) dp_ = g.clone();  
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Solve adjoint equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param_,z,tol);
  Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
  // Build gradient
  ObjectiveBase<Real,std::vector<Real>>::getObjective()->gradient_2(*dp_,state,z,tol);
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(g,adjoint,state,z,tol);
  g.plus(*dp_);
}

template<typename Real>
void C_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (dp_ == nullPtr) dp_ = hv.clone();  
  // Solve state equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param_,z,tol);
  Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
  // Solve adjoint equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param_,z,tol);
  Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
  // Solve state sensitivity equation
  ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity(v,z,tol);
  Vector<Real> &stateSens = *ObjectiveBase<Real,std::vector<Real>>::getStateSens();
  // Solve adjoint sensitivity equation
  ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_sensitivity(v,z,tol);
  Vector<Real> &adjointSens = *ObjectiveBase<Real,std::vector<Real>>::getAdjointSens();
  // Build hessVec
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(hv,adjointSens,state,z,tol);
  ObjectiveBase<Real,std::vector<Real>>::getObjective()->hessVec_21(*dp_,stateSens,state,z,tol);
  hv.plus(*dp_);
  ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointHessian_12(*dp_,adjoint,stateSens,state,z,tol);
  hv.plus(*dp_);
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
