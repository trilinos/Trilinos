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


#ifndef ROL_OED_HOM_OBJECTIVE_BASE_DEF_HPP
#define ROL_OED_HOM_OBJECTIVE_BASE_DEF_HPP

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setConstraint(const Ptr<BilinearConstraint<Real>> &con) {
  con_ = con;
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setObjective(const Ptr<LinearObjective<Real>> &obj) {
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
  state_      = state->clone();
  state_sens_ = state->clone();
  res_        = state->dual().clone();
}

template<typename Real, typename Key>
const Ptr<BilinearConstraint<Real>> ObjectiveBase<Real,Key>::getConstraint() const {
  return con_;
}

template<typename Real, typename Key>
const Ptr<LinearObjective<Real>> ObjectiveBase<Real,Key>::getObjective() const {
  return obj_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getState() const {
  return state_;
}

template<typename Real, typename Key>
const Ptr<Vector<Real>> ObjectiveBase<Real,Key>::getStateSens() const {
  return state_sens_;
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
    // Update equality constraint with new Sim variable.
    if (doUpdate_) con_->update_1(*state_,updateType_,updateIter_);
    // Update full objective function.
    if (obj_ != nullPtr && doUpdate_) obj_->update(*state_,z,updateType_,updateIter_);
    // Store state.
    if (storage_) stateStore_->set(*state_,param);
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
ObjectiveBase<Real,Key>::ObjectiveBase()
  : stateStore_(makePtr<VectorController<Real,Key>>()),
    con_(nullPtr), obj_(nullPtr),
    updateType_(UpdateType::Initial), updateIter_(0),
    storage_(false), doUpdate_(true) {}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::update( const Vector<Real> &z,
                                      UpdateType type,
                                      int iter ) {
  updateType_ = type;
  updateIter_ = iter;
  stateStore_->objectiveUpdate(type);
}

template<typename Real, typename Key>
void ObjectiveBase<Real,Key>::setParameter( const std::vector<Real> &param ) {
  Objective<Real>::setParameter(param);
  if (con_ != nullPtr) con_->setParameter(param);
  if (obj_ != nullPtr) obj_->setParameter(param);
}

} // End Hom Namespace
} // End OED Namespace
} // End ROL Namespace

#endif
