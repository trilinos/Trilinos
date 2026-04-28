// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVED_DEF_HPP
#define ROL_OED_HETOBJECTIVED_DEF_HPP

namespace ROL::OED::Het {

template<typename Real>
void ObjectiveD<Real>::solveState(Vector<Real>& u0, Vector<Real>& u1, Vector<Real>& e, const Vector<Real>& z, unsigned i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? state0Store_->get(u0,i) : false;
  state1Store_->get(u1,i);
  ts_->get(e,{static_cast<Real>(i)});
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M0_->applyInverse(u0,e,z);
    // Solve state equation.
    M1_->applyInverse(u1,e,z);
    // Store state.
    if (storage_) {
      state0Store_->set(u0,i);
      state1Store_->set(u1,i);
    }
  }
}

template<typename Real>
ObjectiveD<Real>::ObjectiveD(const Ptr<MomentOperator<Real>>& M0,
                             const Ptr<MomentOperator<Real>>& M1,
                             const Ptr<const Vector<Real>>& theta,
                             bool storage)
  : M0_(M0), M1_(M1), ts_(makePtr<TraceSampler<Real>>(theta)),
    state0Store_(makePtr<VectorController<Real,unsigned>>()),
    state1Store_(makePtr<VectorController<Real,unsigned>>()),
    u0_(theta->clone()), u1_(u0_->clone()), s_(u0_->clone()),
    e_(theta->dual().clone()), r_(e_->clone()), storage_(storage),
    dim_(theta->dimension()), isDetComputed_(false) {}

template<typename Real>
void ObjectiveD<Real>::update( const Vector<Real>& z, UpdateType type, int iter) {
  M0_->update(z,type,iter);
  M1_->update(z,type,iter);
  state0Store_->objectiveUpdate(type);
  state1Store_->objectiveUpdate(type);
  isDetComputed_ = (type!=UpdateType::Accept ? false : isDetComputed_);
}

template<typename Real>
Real ObjectiveD<Real>::value( const Vector<Real> &z, Real &tol ) {
  if (!isDetComputed_) {
    logdet0_ = M0_->logDeterminant(z);
    logdet1_ = M1_->logDeterminant(z);
    isDetComputed_ = true;
  }
  return logdet0_ - static_cast<Real>(2)*logdet1_;
}

template<typename Real>
void ObjectiveD<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();  
  g.zero();
  for (unsigned i = 0u; i < dim_; ++i) {
    // Solve state equation
    solveState(*u0_,*u1_,*e_,z,i);
    M0_->applySampleMatrices(*g_,*e_,*u0_);
    g.plus(*g_);
    M1_->applySampleMatrices(*g_,*e_,*u1_);
    g.axpy(static_cast<Real>(-2),*g_);
  }
}

template<typename Real>
void ObjectiveD<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = hv.clone();  
  hv.zero();
  for (unsigned i = 0u; i < dim_; ++i) {
    // Solve state equation
    solveState(*u0_,*u1_,*e_,z,i);
    // Solve state sensitivity equation
    M0_->applyDeriv(*r_,*u0_,v);
    M0_->applyInverse(*s_,*r_,z);
    // Assemble Hessian application
    M0_->applySampleMatrices(*g_,*s_,*e_);
    hv.axpy(static_cast<Real>(-1),*g_);
    // Solve state sensitivity equation
    M1_->applyDeriv(*r_,*u1_,v);
    M1_->applyInverse(*s_,*r_,z);
    // Assemble Hessian application
    M1_->applySampleMatrices(*g_,*s_,*e_);
    hv.axpy(static_cast<Real>(2),*g_);
  }
}

} // End ROL::OED::Het Namespace

#endif
