// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEPV_DEF_HPP
#define ROL_OED_HOMOBJECTIVEPV_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
void ObjectivePV<Real>::solveState(PartitionedVector<Real>& u, PartitionedVector<Real>& r, const Vector<Real> &z, const std::vector<Real>& param) {
  // Solve state equation if not done already.
  bool isRhsComputed = storage_ ? rhsStore_->get(r,param) : false;
  bool isStateComputed = storage_ ? stateStore_->get(u,param) : false;
  if (!isStateComputed || !storage_) {
    // Compute rhs if not done already.
    if (!isRhsComputed) {
      for (unsigned i=0; i<nobs_; ++i) {
        ts_->get(*e_,{static_cast<Real>(i)});
        F_->applyAdjoint(*(r.get(i)),*e_,param);
      }
      if (storage_) rhsStore_->set(r,param);
    }
    // Solve and store state equation.
    for (unsigned i=0; i<nobs_; ++i)
      M_->applyInverse(*(u.get(i)),*(r.get(i)),z);
    if (storage_) stateStore_->set(u,param);
  }
}

template<typename Real>
ObjectivePV<Real>::ObjectivePV( const Ptr<MomentOperator<Real>>& M,
                                const Ptr<Factors<Real>>& F,
                                bool storage)
  : M_(M), F_(F), ts_(makePtr<TraceSampler<Real>>(F->createObservationVector(false))),
    stateStore_(makePtr<VectorController<Real,std::vector<Real>>>()),
    rhsStore_(makePtr<VectorController<Real,std::vector<Real>>>()),
    s_(F->createParameterVector(false)), e_(F->createObservationVector(false)),
    d_(F->createParameterVector(true)),
    storage_(storage), nobs_(F->numObservations()) {
  std::vector<Ptr<Vector<Real>>> uvecs(nobs_);
  std::vector<Ptr<Vector<Real>>> rvecs(nobs_);
  for (unsigned i=0; i<nobs_; ++i) {
    uvecs[i] = F->createParameterVector(false);
    rvecs[i] = F->createParameterVector(true);
  }
  u_ = makePtr<PartitionedVector<Real>>(uvecs);
  r_ = makePtr<PartitionedVector<Real>>(rvecs);
}

template<typename Real>
void ObjectivePV<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  M_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectivePV<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  solveState(*u_,*r_,z,Objective<Real>::getParameter());
  // Assemble value
  Real val(0);
  for (unsigned i=0; i<nobs_; ++i)
    val += (r_->get(i))->apply(*(u_->get(i)));
  return val;
}

template<typename Real>
void ObjectivePV<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = g.clone();
  // Solve state equation
  solveState(*u_,*r_,z,Objective<Real>::getParameter());
  // Build gradient
  g.zero();
  for (unsigned i=0; i<nobs_; ++i) {
    M_->applySampleMatrices(*g_,*(u_->get(i)),*(u_->get(i)));
    g.plus(*g_);
  }
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void ObjectivePV<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = hv.clone();
  // Solve state equation
  solveState(*u_,*r_,z,Objective<Real>::getParameter());
  hv.zero();
  // Assemble hessVec
  for (unsigned i=0; i<nobs_; ++i) {
    // Solve state sensitivity equation
    M_->applyDeriv(*d_,*(u_->get(i)),v);
    M_->applyInverse(*s_,*d_,z);
    M_->applySampleMatrices(*g_,*s_,*(u_->get(i)));
    hv.plus(*g_);
  }
  hv.scale(static_cast<Real>(2));
}

} // END ROL::OED::Hom Namespace

#endif
