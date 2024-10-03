// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_A_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_A_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
A_Objective<Real>::A_Objective( const Ptr<BilinearConstraint<Real>>   &con,
             const Ptr<QuadraticObjective<Real>> &obj,
             const Ptr<Vector<Real>>             &state,
             const std::vector<Real>             &weight,
             bool storage) : weight_(weight) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real A_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  const int dim = weight_.size();
  Real val(0);
  std::vector<Real> param(1);
  for (int i = 0; i < dim; ++i) {
    param[0] = static_cast<Real>(i);
    // Solve state equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param,z,tol);
    Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
    // Get objective function value
    val += weight_[i]*ObjectiveBase<Real,std::vector<Real>>::getObjective()->value(state,z,tol);
  }
  return val;
}

template<typename Real>
void A_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();
  const int dim = weight_.size();
  g.zero();
  std::vector<Real> param(1);
  for (int i = 0; i < dim; ++i) {
    param[0] = static_cast<Real>(i);
    // Solve state equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param,z,tol);
    Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
    // Solve adjoint equation
    ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param,z,tol);
    Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
    // Build gradient
    ObjectiveBase<Real,std::vector<Real>>::getObjective()->gradient_2(*g_,state,z,tol);
    g.axpy(weight_[i],*g_);
    ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(*g_,adjoint,state,z,tol);
    g.axpy(weight_[i],*g_);
  }
}

template<typename Real>
void A_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = hv.clone();
  const int dim = weight_.size();
  hv.zero();
  std::vector<Real> param(1);
  for (int i = 0; i < dim; ++i) {
    param[0] = static_cast<Real>(i);
    // Solve state equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param,z,tol);
    Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
    // Solve adjoint equation
    ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation(param,z,tol);
    Vector<Real> &adjoint = *ObjectiveBase<Real,std::vector<Real>>::getAdjoint();
    // Solve state sensitivity equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity(v,z,tol);
    Vector<Real> &stateSens = *ObjectiveBase<Real,std::vector<Real>>::getStateSens();
    // Solve adjoint sensitivity equation
    ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_sensitivity(v,z,tol);
    Vector<Real> &adjointSens = *ObjectiveBase<Real,std::vector<Real>>::getAdjointSens();
    // Build hessVec
    ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(*g_,adjointSens,state,z,tol);
    hv.axpy(weight_[i],*g_);
    ObjectiveBase<Real,std::vector<Real>>::getObjective()->hessVec_21(*g_,stateSens,state,z,tol);
    hv.axpy(weight_[i],*g_);
    ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointHessian_12(*g_,adjoint,stateSens,state,z,tol);
    hv.axpy(weight_[i],*g_);
  }
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
