// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_REDUCED_OBJECTIVE_SIMOPT_DEF_H
#define ROL_REDUCED_OBJECTIVE_SIMOPT_DEF_H

namespace ROL {

template<typename Real>
Reduced_Objective_SimOpt<Real>::Reduced_Objective_SimOpt(
    const Ptr<Objective_SimOpt<Real>> &obj,
    const Ptr<Constraint_SimOpt<Real>> &con,
    const Ptr<Vector<Real>> &state,
    const Ptr<Vector<Real>> &control,
    const Ptr<Vector<Real>> &adjoint,
    const bool storage,
    const bool useFDhessVec)
  : obj_(obj), con_(con),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), ngrad_(0), nhess_(0), nprec_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false) {
  stateStore_   = makePtr<VectorController<Real>>();
  adjointStore_ = makePtr<VectorController<Real>>();
  state_        = state->clone(); state_->set(*state);
  adjoint_      = adjoint->clone();
  state_sens_   = state->clone();
  adjoint_sens_ = adjoint->clone();
  dualstate_    = state->dual().clone();
  dualstate1_   = state->dual().clone();
  dualadjoint_  = adjoint->dual().clone();
  dualcontrol_  = control->dual().clone();
}

template<typename Real>
Reduced_Objective_SimOpt<Real>::Reduced_Objective_SimOpt(
    const Ptr<Objective_SimOpt<Real>> &obj,
    const Ptr<Constraint_SimOpt<Real>> &con,
    const Ptr<Vector<Real>> &state,
    const Ptr<Vector<Real>> &control,
    const Ptr<Vector<Real>> &adjoint,
    const Ptr<Vector<Real>> &dualstate,
    const Ptr<Vector<Real>> &dualcontrol,
    const Ptr<Vector<Real>> &dualadjoint,
    const bool storage,
    const bool useFDhessVec)
  : obj_(obj), con_(con),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), ngrad_(0), nhess_(0), nprec_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false) {
  stateStore_   = makePtr<VectorController<Real>>();
  adjointStore_ = makePtr<VectorController<Real>>();
  state_        = state->clone(); state_->set(*state);
  adjoint_      = adjoint->clone();
  state_sens_   = state->clone();
  adjoint_sens_ = adjoint->clone();
  dualstate_    = dualstate->clone();
  dualstate1_   = dualstate->clone();
  dualadjoint_  = dualadjoint->clone();
  dualcontrol_  = dualcontrol->clone();
}

template<typename Real>
Reduced_Objective_SimOpt<Real>::Reduced_Objective_SimOpt(
    const Ptr<Objective_SimOpt<Real>> &obj,
    const Ptr<Constraint_SimOpt<Real>> &con,
    const Ptr<VectorController<Real>> &stateStore,
    const Ptr<Vector<Real>> &state,
    const Ptr<Vector<Real>> &control,
    const Ptr<Vector<Real>> &adjoint,
    const bool storage,
    const bool useFDhessVec)
  : obj_(obj), con_(con), stateStore_(stateStore),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), ngrad_(0), nhess_(0), nprec_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false) {
  adjointStore_ = makePtr<VectorController<Real>>();
  state_        = state->clone(); state_->set(*state);
  adjoint_      = adjoint->clone();
  state_sens_   = state->clone();
  adjoint_sens_ = adjoint->clone();
  dualstate_    = state->dual().clone();
  dualstate1_   = state->dual().clone();
  dualadjoint_  = adjoint->dual().clone();
  dualcontrol_  = control->dual().clone();
}

template<typename Real>
Reduced_Objective_SimOpt<Real>::Reduced_Objective_SimOpt(
    const Ptr<Objective_SimOpt<Real>> &obj,
    const Ptr<Constraint_SimOpt<Real>> &con,
    const Ptr<VectorController<Real>> &stateStore,
    const Ptr<Vector<Real>> &state,
    const Ptr<Vector<Real>> &control,
    const Ptr<Vector<Real>> &adjoint,
    const Ptr<Vector<Real>> &dualstate,
    const Ptr<Vector<Real>> &dualcontrol,
    const Ptr<Vector<Real>> &dualadjoint,
    const bool storage,
    const bool useFDhessVec)
  : obj_(obj), con_(con), stateStore_(stateStore),
    storage_(storage), useFDhessVec_(useFDhessVec),
    nupda_(0), nvalu_(0), ngrad_(0), nhess_(0), nprec_(0),
    nstat_(0), nadjo_(0), nssen_(0), nasen_(0),
    updateFlag_(true), updateIter_(0), updateType_(UpdateType::Initial),
    newUpdate_(false) {
  adjointStore_ = makePtr<VectorController<Real>>();
  state_        = state->clone(); state_->set(*state);
  adjoint_      = adjoint->clone();
  state_sens_   = state->clone();
  adjoint_sens_ = adjoint->clone();
  dualstate_    = dualstate->clone();
  dualstate1_   = dualstate->clone();
  dualadjoint_  = dualadjoint->clone();
  dualcontrol_  = dualcontrol->clone();
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::update( const Vector<Real> &z, bool flag, int iter ) {
  nupda_++;
  newUpdate_  = false;
  updateFlag_ = flag;
  updateIter_ = iter;
  stateStore_->objectiveUpdate(true);
  adjointStore_->objectiveUpdate(flag);
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::update( const Vector<Real> &z, UpdateType type, int iter ) {
  nupda_++;
  newUpdate_  = true;
  updateType_ = type;
  updateIter_ = iter;
  stateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real>
Real Reduced_Objective_SimOpt<Real>::value( const Vector<Real> &z, Real &tol ) {
  nvalu_++;
  // Solve state equation
  solve_state_equation(z,tol);
  // Get objective function value
  return obj_->value(*state_,z,tol);
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  ngrad_++;
  // Solve state equation
  solve_state_equation(z,tol);
  // Solve adjoint equation
  solve_adjoint_equation(z,tol);
  // Evaluate the full gradient wrt z
  obj_->gradient_2(*dualcontrol_,*state_,z,tol);
  // Build gradient
  con_->applyAdjointJacobian_2(g,*adjoint_,*state_,z,tol);
  g.plus(*dualcontrol_);
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  nhess_++;
  if ( useFDhessVec_ ) {
    Objective<Real>::hessVec(hv,v,z,tol);
  }
  else {
    // Solve state equation
    solve_state_equation(z,tol);
    // Solve adjoint equation
    solve_adjoint_equation(z,tol);
    // Solve state sensitivity equation
    solve_state_sensitivity(v,z,tol);
    // Solve adjoint sensitivity equation
    solve_adjoint_sensitivity(v,z,tol);
    // Build hessVec
    con_->applyAdjointJacobian_2(hv,*adjoint_sens_,*state_,z,tol);
    obj_->hessVec_21(*dualcontrol_,*state_sens_,*state_,z,tol);
    hv.plus(*dualcontrol_);
    obj_->hessVec_22(*dualcontrol_,v,*state_,z,tol);
    hv.plus(*dualcontrol_);
    con_->applyAdjointHessian_12(*dualcontrol_,*adjoint_,*state_sens_,*state_,z,tol);
    hv.plus(*dualcontrol_);
    con_->applyAdjointHessian_22(*dualcontrol_,*adjoint_,v,*state_,z,tol);
    hv.plus(*dualcontrol_);
  }
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  nprec_++;
  Pv.set(v.dual());
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  con_->setParameter(param);
  obj_->setParameter(param);
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::summarize(std::ostream &stream) const {
  stream << std::endl;
  stream << std::string(80,'=') << std::endl;
  stream << "  ROL::Reduced_Objective_SimOpt::summarize" << std::endl;
  stream << "    Number of calls to update:            " << nupda_ << std::endl;
  stream << "    Number of calls to value:             " << nvalu_ << std::endl;
  stream << "    Number of calls to gradient:          " << ngrad_ << std::endl;
  stream << "    Number of calls to hessvec:           " << nhess_ << std::endl;
  stream << "    Number of calls to precond:           " << nprec_ << std::endl;
  stream << "    Number of state solves:               " << nstat_ << std::endl;
  stream << "    Number of adjoint solves:             " << nadjo_ << std::endl;
  stream << "    Number of state sensitivity solves:   " << nssen_ << std::endl;
  stream << "    Number of adjoint sensitivity solves: " << nssen_ << std::endl;
  stream << std::string(80,'=') << std::endl;
  stream << std::endl;
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::reset() {
  nupda_ = 0; nvalu_ = 0; ngrad_ = 0; nhess_ = 0; nprec_ = 0;
  nstat_ = 0; nadjo_ = 0; nssen_ = 0; nasen_ = 0;
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::solve_state_equation(const Vector<Real> &z, Real &tol) {
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(*state_,Objective<Real>::getParameter()) : false;
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Update equality constraint with new Opt variable.
    if (newUpdate_) con_->update_2(z,updateType_,updateIter_);
    else            con_->update_2(z,updateFlag_,updateIter_);
    // Solve state equation.
    con_->solve(*dualadjoint_,*state_,z,tol);
    nstat_++;
    // Update equality constraint with new Sim variable.
    if (newUpdate_) con_->update_1(*state_,updateType_,updateIter_);
    else            con_->update_1(*state_,updateFlag_,updateIter_);
    // Update full objective function.
    if (newUpdate_) obj_->update(*state_,z,updateType_,updateIter_);
    else            obj_->update(*state_,z,updateFlag_,updateIter_);
    // Store state.
    if (storage_)   stateStore_->set(*state_,Objective<Real>::getParameter());
  }
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::solve_adjoint_equation(const Vector<Real> &z, Real &tol) {
  // Check if adjoint has been computed.
  bool isComputed = storage_ ? adjointStore_->get(*adjoint_,Objective<Real>::getParameter()) : false;
  // Solve adjoint equation if not done already.
  if (!isComputed || !storage_) {
    // Evaluate the full gradient wrt u
    obj_->gradient_1(*dualstate_,*state_,z,tol);
    // Solve adjoint equation
    con_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,z,tol);
    adjoint_->scale(static_cast<Real>(-1));
    nadjo_++;
    // Store adjoint
    if (storage_) adjointStore_->set(*adjoint_,Objective<Real>::getParameter());
  }
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  // Solve state sensitivity equation
  con_->applyJacobian_2(*dualadjoint_,v,*state_,z,tol);
  dualadjoint_->scale(static_cast<Real>(-1));
  con_->applyInverseJacobian_1(*state_sens_,*dualadjoint_,*state_,z,tol);
  nssen_++;
}

template<typename Real>
void Reduced_Objective_SimOpt<Real>::solve_adjoint_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  // Evaluate full hessVec in the direction (s,v)
  obj_->hessVec_11(*dualstate_,*state_sens_,*state_,z,tol);
  obj_->hessVec_12(*dualstate1_,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Apply adjoint Hessian of constraint
  con_->applyAdjointHessian_11(*dualstate1_,*adjoint_,*state_sens_,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  con_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Solve adjoint sensitivity equation
  dualstate_->scale(static_cast<Real>(-1));
  con_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,z,tol);
  nasen_++;
}

} // namespace ROL

#endif
