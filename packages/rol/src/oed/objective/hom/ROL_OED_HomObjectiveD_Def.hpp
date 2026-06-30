// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVED_DEF_HPP
#define ROL_OED_HOMOBJECTIVED_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
void ObjectiveD<Real>::solveState(Vector<Real>& u, Vector<Real>& e, const Vector<Real> &z, int i) {
  std::vector<Real> param = {static_cast<Real>(i)};
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(u,i) : false;
  ts_->get(e,param);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M_->applyInverse(u,e,z);
    // Store state.
    if (storage_) stateStore_->set(u,i);
  }
}

template<typename Real>
ObjectiveD<Real>::ObjectiveD( const Ptr<MomentOperator<Real>>& M,
                              const Ptr<const Vector<Real>>& theta,
                              bool storage)
  : M_(M), ts_(makePtr<TraceSampler<Real>>(theta)),
    stateStore_(makePtr<VectorController<Real,int>>()),
    u_(theta->clone()), s_(u_->clone()), e_(theta->dual().clone()), r_(e_->clone()),
    storage_(storage), dim_(theta->dimension()), isDetComputed_(false) {}

template<typename Real>
void ObjectiveD<Real>::update(const Vector<Real> &z, UpdateType type, int iter) {
  M_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
  isDetComputed_ = (type!=UpdateType::Accept ? false : isDetComputed_);
}

template<typename Real>
Real ObjectiveD<Real>::value( const Vector<Real> &z, Real &tol ) {
  if (!isDetComputed_) {
    logdet_ = M_->logDeterminant(z);
    isDetComputed_ = true;
  }
  return -logdet_;
}

template<typename Real>
void ObjectiveD<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();
  g.zero();
  for (int i = 0; i < dim_; ++i) {
    // Solve state equation
    solveState(*u_,*e_,z,i);
    M_->applySampleMatrices(*g_,*e_,*u_);
    g.plus(*g_);
  }
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void ObjectiveD<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = hv.clone();
  hv.zero();
  for (int i = 0; i < dim_; ++i) {
    // Solve state equation
    solveState(*u_,*e_,z,i);
    // Solve state sensitivity equation
    M_->applyDeriv(*r_,*u_,v);
    M_->applyInverse(*s_,*r_,z);
    // Assemble Hessian application
    M_->applySampleMatrices(*g_,*s_,*e_);
    hv.plus(*g_);
  }
}

} // END ROL::OED::Hom Namespace

#endif
