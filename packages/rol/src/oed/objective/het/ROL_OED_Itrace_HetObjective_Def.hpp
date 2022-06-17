// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package Copyright (2014)
//               Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive license
// for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the contributors may
// be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL SANDIA CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers: Drew Kouri   (dpkouri@sandia.gov) and
// Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_OED_ITRACE_HET_OBJECTIVE_DEF_HPP
#define ROL_OED_ITRACE_HET_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
void Itrace_Objective<Real>::solveStateEquation(Vector<Real> &state, const Vector<Real> &u,
                        const Vector<Real> &z, int i, Real &tol) {
  bool isComputed = (storage_ ? stateStore_->get(state,i) : false);
  if (!isComputed) {
    getConstraint()->applyInverseJacobian_1(state,*b_[i],u,z,tol);
    if (storage_) stateStore_->set(state,i);
  }
}

template<typename Real>
void Itrace_Objective<Real>::solveAdjointEquation(Vector<Real> &adjoint, const Vector<Real> &state,
                          const Vector<Real> &u, const Vector<Real> &z,
                          VectorController<Real,int> &store,
                          int i, Real &tol) {
  bool isComputed = (storage_ ? store.get(adjoint,i) : false);
  if (!isComputed) {
    M_->applyJacobian_1(*rhs_,state,u,z,tol);
    getConstraint()->applyInverseJacobian_1(adjoint,*rhs_,u,z,tol);
    if (storage_) store.set(adjoint,i);
  }
}

template<typename Real>
void Itrace_Objective<Real>::solveStateSensitivityEquation(Vector<Real> &sens,
                             const Vector<Real> &v, const Vector<Real> &u,
                             const Vector<Real> &z, Real &tol) {
  getConstraint()->applyJacobian_2(*rhs_,v,u,z,tol);
  getConstraint()->applyInverseJacobian_1(sens,*rhs_,u,z,tol);
}

template<typename Real>
void Itrace_Objective<Real>::solveAdjointSensitivityEquation(Vector<Real> &sens,
                               const Vector<Real> &v,
                               const Vector<Real> &p, const Vector<Real> &s,
                               const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
  M_->applyJacobian_1(*rhs_,s,u,z,tol);
  M_->applyJacobian_2(*rhs1_,v,u,z,tol);
  rhs_->axpy(static_cast<Real>(-1),*rhs1_);
  getConstraint()->applyJacobian_2(*rhs1_,v,p,z,tol);
  rhs_->plus(*rhs1_);
  getConstraint()->applyInverseJacobian_1(sens,*rhs_,u,z,tol);
}

template<typename Real>
Itrace_Objective<Real>::Itrace_Objective( const Ptr<BilinearConstraint<Real>>     &con,
                  const Ptr<QuadraticObjective<Real>>   &obj,
                  const Ptr<Vector<Real>>               &state,
                  const Ptr<SampleGenerator<Real>>      &sampler,
                  const std::vector<Real>               &weight,
                  bool storage)
  : M_(obj->getM()), sampler_(sampler), weight_(weight),
    state_(state->clone()), adjoint1_(state->clone()),
    adjoint2_(state->clone()), sens1_(state->clone()),
    sens2_(state->clone()), sadj1_(state->clone()),
    sadj2_(state->clone()), rhs_(state->dual().clone()),
    rhs1_(state->dual().clone()),
    stateStore_(makePtr<VectorController<Real,int>>()),
    adjoint1Store_(makePtr<VectorController<Real,int>>()),
    adjoint2Store_(makePtr<VectorController<Real,int>>()),
    storage_(storage) {
  ObjectiveBase<Real,std::vector<Real>>::setConstraint(con);
  ObjectiveBase<Real,std::vector<Real>>::setObjective(obj);
  ObjectiveBase<Real,std::vector<Real>>::setStorage(storage);
  ObjectiveBase<Real,std::vector<Real>>::initialize(state);

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
             int iter) {
  ObjectiveBase<Real,std::vector<Real>>::update(z,type,iter);
  stateStore_->objectiveUpdate(type);
  adjoint1Store_->objectiveUpdate(type);
  adjoint2Store_->objectiveUpdate(type);
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
    solveStateEquation(*state_,state,z,i,tol);
    // Get objective function value
    M_->applyJacobian_1(*rhs_,*state_,state,z,tol);
    val += weight_[i]*rhs_->apply(state);
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
    solveStateEquation(*state_,state,z,i,tol);
    // Solve adjoint equation
    solveAdjointEquation(*adjoint1_,*state_,state,z,*adjoint1Store_,i,tol);
    solveAdjointEquation(*adjoint2_, state, state,z,*adjoint2Store_,i,tol);
    // Build gradient
    M_->applyAdjointJacobian_2(*Xa_,*state_,state,z,tol);
    g.axpy(weight_[i],*Xa_);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*adjoint1_,state,z,tol);
    g.axpy(-weight_[i],*Xa_);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*adjoint2_,*state_,z,tol);
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
    solveStateEquation(*state_,state,z,i,tol);
    // Solve adjoint equation
    solveAdjointEquation(*adjoint1_,*state_,state,z,*adjoint1Store_,i,tol);
    solveAdjointEquation(*adjoint2_, state, state,z,*adjoint2Store_,i,tol);
    // Solve state sensitivity equation
    solveStateSensitivityEquation(*sens1_,v,state,z,tol);
    solveStateSensitivityEquation(*sens2_,v,*state_,z,tol);
    // Solve adjoint sensitivity equation
    solveAdjointSensitivityEquation(*sadj1_,v,*adjoint1_,*sens2_,*state_,z,tol);
    solveAdjointSensitivityEquation(*sadj2_,v,*adjoint2_,*sens1_,state,z,tol);
    // Build hessvec
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sadj1_,state,z,tol);
    hv.axpy(weight_[i],*Xa_);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sadj2_,*state_,z,tol);
    hv.axpy(weight_[i],*Xa_);
    M_->applyAdjointJacobian_2(*Xa_,*sens1_,*state_,z,tol);
    hv.axpy(-weight_[i],*Xa_);
    M_->applyAdjointJacobian_2(*Xa_,*sens2_,state,z,tol);
    hv.axpy(-weight_[i],*Xa_);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sens1_,*adjoint1_,z,tol);
    hv.axpy(weight_[i],*Xa_);
    getConstraint()->applyAdjointJacobian_2(*Xa_,*sens2_,*adjoint2_,z,tol);
    hv.axpy(weight_[i],*Xa_);
  }
}

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#endif
