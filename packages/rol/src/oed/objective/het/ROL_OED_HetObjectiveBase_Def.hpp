// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HET_OBJECTIVE_BASE_DEF_HPP
#define ROL_OED_HET_OBJECTIVE_BASE_DEF_HPP

namespace ROL {
namespace OED {
namespace Het {

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setConstraint(const Ptr<BilinearConstraint<Real>> &con) {
  con_ = con;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setObjective(const Ptr<QuadraticObjective<Real>> &obj) {
  obj_ = obj;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setStorage(bool storage) {
  storage_ = storage;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setUpdate(bool doUpdate) {
  doUpdate_ = doUpdate;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::initialize(const Ptr<Vector<Real>> &state) {
  state_        = state->clone();
  state_sens_   = state->clone();
  adjoint_      = state->clone();
  adjoint_sens_ = state->clone();
  dualstate_    = state->dual().clone();
  dualstate1_   = state->dual().clone();
  res_          = state->dual().clone();
}

template<typename Real, typename Key>
const Ptr<BilinearConstraint<Real>> ObjectiveBase<Real,Key>::getConstraint() const {
  return con_;
}

template<typename Real, typename Key>
const Ptr<QuadraticObjective<Real>> ObjectiveBase<Real,Key>::getObjective() const {
  return obj_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getState() const {
  return state_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getAdjoint() const {
  return adjoint_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getStateSens() const {
  return state_sens_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getAdjointSens() const {
  return adjoint_sens_;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::solve_state_equation(const Key &param,
                          const Vector<Real> &z,
                                Real &tol) { 
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(*state_,param) : false;
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    if (obj_ != nullPtr) obj_->setParameter(param);
    con_->setParameter(param);
    // Update equality constraint with new Opt variable.
    if (doUpdate_) con_->update_2(z,updateType_,updateIter_);
    // Solve state equation.
    con_->solve(*res_,*state_,z,tol);
    //// Update equality constraint with new Sim variable.
    if (doUpdate_) con_->update_1(*state_,updateType_,updateIter_);
    // Update full objective function.
    if (obj_ != nullPtr && doUpdate_) obj_->update(*state_,z,updateType_,updateIter_);
    // Store state.
    if (storage_) stateStore_->set(*state_,param);
  }
  if (doUpdate_) con_->update_1(*state_,updateType_,updateIter_);
  if (obj_ != nullPtr && doUpdate_) obj_->update(*state_,z,updateType_,updateIter_);
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::solve_adjoint_equation(const Key &param,
                            const Vector<Real> &z,
                                  Real &tol) { 
  // Check if adjoint has been computed.
  bool isComputed = false;
  if (storage_) {
    isComputed = adjointStore_->get(*adjoint_,param);
  }
  // Solve adjoint equation if not done already.
  if (!isComputed || !storage_) {
    obj_->setParameter(param);
    con_->setParameter(param);
    // Evaluate the full gradient wrt u
    obj_->gradient_1(*dualstate_,*state_,z,tol);
    // Solve adjoint equation
    con_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,z,tol);
    adjoint_->scale(static_cast<Real>(-1));
    // Store adjoint
    if (storage_) {
      adjointStore_->set(*adjoint_,param);
    }
  }
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::solve_state_sensitivity(const Vector<Real> &v,
                             const Vector<Real> &z,
                                   Real &tol) {
  // Solve state sensitivity equation
  con_->applyJacobian_2(*res_,v,*state_,z,tol);
  res_->scale(static_cast<Real>(-1));
  con_->applyInverseJacobian_1(*state_sens_,*res_,*state_,z,tol);
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::solve_adjoint_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
  // Evaluate full hessVec in the direction (s,v)
  obj_->hessVec_11(*dualstate_,*state_sens_,*state_,z,tol);
  obj_->hessVec_12(*dualstate1_,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Apply adjoint Hessian of constraint
  con_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,z,tol);
  dualstate_->plus(*dualstate1_);
  // Solve adjoint sensitivity equation
  dualstate_->scale(static_cast<Real>(-1));
  con_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,z,tol);
}

template<typename Real, typename Key>
ObjectiveBase<Real,Key>::ObjectiveBase()
  : con_(nullPtr), obj_(nullPtr),
    storage_(false), doUpdate_(true), updateType_(UpdateType::Initial), updateIter_(0) {
  stateStore_   = makePtr<VectorController<Real,Key>>();
  adjointStore_ = makePtr<VectorController<Real,Key>>();
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::update( const Vector<Real> &z,
             UpdateType type,
             int iter ) {
  updateType_ = type;
  updateIter_ = iter;
  stateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setParameter( const std::vector<Real> &param ) {
  Objective<Real>::setParameter(param);
  con_->setParameter(param);
  obj_->setParameter(param);
}

} // End Het Namespace
} // End OED Namespace
} // End ROL Namespace

#endif
