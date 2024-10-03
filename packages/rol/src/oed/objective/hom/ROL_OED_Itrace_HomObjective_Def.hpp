// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_ITRACE_HOM_OBJECTIVE_DEF_HPP
#define ROL_OED_ITRACE_HOM_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
void Itrace_Objective<Real>::solveAdjointEquation(Vector<Real> &adjoint, const Vector<Real> &u,
                          const Vector<Real> &z, int i, Real &tol) {
  bool isComputed = (storage_ ? adjointStore_->get(adjoint,i) : false);
  if (!isComputed) {
    getConstraint()->applyInverseJacobian_1(adjoint,*b_[i],u,z,tol);
    if (storage_) adjointStore_->set(adjoint,i);
  }
}

template<typename Real>
Itrace_Objective<Real>::Itrace_Objective( const Ptr<BilinearConstraint<Real>> &con,
                  const Ptr<LinearObjective<Real>>  &obj,
                  const Ptr<Vector<Real>>           &state,
                  const Ptr<SampleGenerator<Real>>  &sampler,
                  const std::vector<Real>           &weight,
                  bool storage)
  : sampler_(sampler), weight_(weight),
    adjoint_(state->clone()), rhs_(state->dual().clone()), sens_(state->clone()),
    adjointStore_(makePtr<VectorController<Real,int>>()),
    storage_(storage) {
  setConstraint(con);
  setObjective(obj);
  setStorage(storage);
  initialize(state);

  int dim = weight_.size();
  Ptr<Vector<Real>> F = state->dual().clone(),
                              g = state->dual().clone(),
                              b = state->dual().clone();
  b_.clear(); b_.resize(dim);
  for (int i = 0; i < dim; ++i) {
    b->zero();
    b_[i] = state->dual().clone();
    getConstraint()->getTraceSample(*g,{static_cast<Real>(i)});
    for (int j = 0; j < sampler_->numMySamples(); ++j) {
      getConstraint()->getFactor(*F,sampler_->getMyPoint(j));
      b->axpy(F->dot(*g)*sampler_->getMyWeight(j),*F);
    }
    sampler_->sumAll(*b,*b_[i]);
  }
}

template<typename Real>
void Itrace_Objective<Real>::update( const Vector<Real> &z,
             UpdateType type,
             int iter ) {
  ObjectiveBase<Real,std::vector<Real>>::update(z,type,iter);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real>
Real Itrace_Objective<Real>::value( const Vector<Real> &z, Real &tol ) {
  const int dim = weight_.size();
  Real val(0);
  for (int i = 0; i < dim; ++i) {
    // Solve state equation
    std::vector<Real> param = {static_cast<Real>(i)};
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    // Get objective function value
    val += weight_[i]*b_[i]->apply(state);
  }
  return val;
}

template<typename Real>
void Itrace_Objective<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (Xa_==nullPtr) Xa_ = z.clone();
  g.zero();
  const int dim = weight_.size();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    // Solve adjoint equation
    solveAdjointEquation(*adjoint_,state,z,i,tol);
    // Build gradient
    getConstraint()->applyAdjointJacobian_2(*Xa_,*adjoint_,state,z,tol);
    g.axpy(-weight_[i],*Xa_);
  }
}

template<typename Real>
void Itrace_Objective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (Xa_==nullPtr) Xa_ = z.clone();
  const int dim = weight_.size();
  hv.zero();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solve_state_equation(param,z,tol);
    Vector<Real> &state = *getState();
    // Solve adjoint equation
    solveAdjointEquation(*adjoint_,state,z,i,tol);
    // Solve state and adjoint sensitivity equation
    getConstraint()->applyJacobian_2(*rhs_,v,state,z,tol);
    getConstraint()->applyInverseJacobian_1(*sens_,*rhs_,state,z,tol);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sens_,*adjoint_,z,tol);
    hv.axpy(weight_[i],*Xa_);
    getConstraint()->applyJacobian_2(*rhs_,v,*adjoint_,z,tol);
    getConstraint()->applyInverseJacobian_1(*sens_,*rhs_,state,z,tol);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sens_,state,z,tol);
    hv.axpy(weight_[i],*Xa_);
  }
}

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
