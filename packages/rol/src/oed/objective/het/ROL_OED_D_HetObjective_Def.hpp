// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_D_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_D_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
void D_Objective<Real>::solve(int key, const Vector<Real> &z, Real &tol) { 
  // Solve state equation if not done already.
  if (!ustore_->get(*u_,key)) {
    con0_->setParameter({static_cast<Real>(key)});
    con0_->solve(*g_,*u_,z,tol);
    ustore_->set(*u_,key); // Store state.
  }
}

template<typename Real>
void D_Objective<Real>::solveSens(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  con0_->applyJacobian_2(*g_,v,*u_,z,tol);
  g_->scale(static_cast<Real>(-1));
  con0_->applyInverseJacobian_1(*u_,*g_,*u_,z,tol);
}

template<typename Real>
D_Objective<Real>::D_Objective( const Ptr<BilinearConstraint<Real>> &con0, // Objective Moment Operator
             const Ptr<BilinearConstraint<Real>> &con1, // Constraint Moment Operator
             const Ptr<Vector<Real>>           &state,
             bool storage)
  : con0_(con0), ustore_(makePtr<VectorController<Real,int>>()),
    u_(state->clone()), g_(state->dual().clone()), dim_(state->dimension()),
    isDetComputed_(false) {
  setConstraint(con1);
  setStorage(storage);
  setUpdate(false);
  initialize(state);
}

template<typename Real>
void D_Objective<Real>::update(const Vector<Real> &z, UpdateType type, int iter) {
  ObjectiveBase<Real,std::vector<Real>>::update(z,type,iter);
  con0_->update_2(z,type,iter);
  getConstraint()->update_2(z,type,iter);
  isDetComputed_ = (type!=UpdateType::Accept ? false : isDetComputed_);
  ustore_->objectiveUpdate(type);
}

template<typename Real>
Real D_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  if (!isDetComputed_) {
    det0_ = con0_->logDeterminant(z);
    det1_ = ObjectiveBase<Real,std::vector<Real>>::getConstraint()->logDeterminant(z);
    isDetComputed_ = true;
  }
  return det0_-static_cast<Real>(2)*det1_;
}

template<typename Real>
void D_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (p_ == nullPtr) p_ = g.clone();
  g.zero();
  for (int i = 0; i < dim_; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    getConstraint()->applyAdjointJacobian_2(*p_,*state.basis(i),state,z,tol);
    g.axpy(static_cast<Real>(-2),*p_);

    // Solve state equation
    solve(i,z,tol);
    con0_->applyAdjointJacobian_2(*p_,*u_->basis(i),*u_,z,tol);
    g.plus(*p_);
  }
}

template<typename Real>
void D_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (p_ == nullPtr) p_ = hv.clone();
  hv.zero();
  for (int i = 0; i < dim_; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    // Solve state sensitivity equation
    solve_state_sensitivity(v,z,tol);
    Vector<Real> &sens = *getStateSens();
    getConstraint()->applyAdjointJacobian_2(*p_,*state.basis(i),sens,z,tol);
    hv.axpy(static_cast<Real>(-2),*p_);

    // Solve state equation
    solve(i,z,tol);
    // Solve state sensitivity equation
    solveSens(v,z,tol);
    con0_->applyAdjointJacobian_2(*p_,*u_->basis(i),*u_,z,tol);
    hv.plus(*p_);
  }
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
