// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_A_HOM_OBJECTIVE_DEF_HPP
#define ROL_OED_A_HOM_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
void A_Objective<Real>::initialize(const Vector<Real> &p) {
  if (!isInitialized_) {
    g_ = p.dual().clone();
    isInitialized_ = true;
  }
}

template<typename Real>
A_Objective<Real>::A_Objective( const Ptr<BilinearConstraint<Real>> &con,
             const Ptr<LinearObjective<Real>>  &obj,
             const Ptr<Vector<Real>>           &state,
             const std::vector<Real>           &weight,
             bool storage) : weight_(weight), isInitialized_(false) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);
}

template<typename Real>
Real A_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  const int dim = weight_.size();
  Real val(0);
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
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
  initialize(z);
  const int dim = weight_.size();
  g.zero();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param,z,tol);
    Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
    // Build gradient
    ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(*g_,state,state,z,tol);
    g.axpy(-weight_[i],*g_);
  }
}

template<typename Real>
void A_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  initialize(z);
  const int dim = weight_.size();
  hv.zero();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_equation(param,z,tol);
    Vector<Real> &state = *ObjectiveBase<Real,std::vector<Real>>::getState();
    // Solve state sensitivity equation
    ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity(v,z,tol);
    Vector<Real> &stateSens = *ObjectiveBase<Real,std::vector<Real>>::getStateSens();
    // Build hessVec
    ObjectiveBase<Real,std::vector<Real>>::getConstraint()->applyAdjointJacobian_2(*g_,stateSens,state,z,tol);
    hv.axpy(static_cast<Real>(-2)*weight_[i],*g_);
  }
}

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
