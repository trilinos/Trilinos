// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEPV_DEF_HPP
#define ROL_OED_HETOBJECTIVEPV_DEF_HPP

namespace ROL::OED::Het {

template<typename Real>
void ObjectivePV<Real>::solveState(PartitionedVector<Real>& u, PartitionedVector<Real>& Mu, const Vector<Real>& z, const std::vector<Real>& param) {
  // Solve state equation if not done already.
  bool isStateComputed = storage_ ? stateStore_->get(u,param) : false;
  MstateStore_->get(Mu,param);
  if (!isStateComputed || !storage_) {
    // Compute rhs if not done already.
    bool isRhsComputed = storage_ ? rhsStore_->get(*r_,param) : false;
    if (!isRhsComputed) {
      for (unsigned i=0u; i<nobs_; ++i) {
        ts_->get(*e_,{static_cast<Real>(i)});
        F_->applyAdjoint(*(r_->get(i)),*e_,param);
      }
      if (storage_) rhsStore_->set(*r_,param);
    }
    // Solve and store state equation.
    for (unsigned i=0u; i<nobs_; ++i) {
      M1_->applyInverse(*(u.get(i)),*(r_->get(i)),z);
      M0_->apply(*(Mu.get(i)),*(u.get(i)),z);
    }
    if (storage_) {
      stateStore_->set(u,param);
      MstateStore_->set(Mu,param);
    }
  }
}

template<typename Real>
void ObjectivePV<Real>::solveAdjoint(PartitionedVector<Real>& p, const PartitionedVector<Real>& Mu, const Vector<Real>& z, const std::vector<Real>& param) {
  // Solve adjoint equation if not done already.
  bool isAdjointComputed = storage_ ? adjointStore_->get(p,param) : false;
  if (!isAdjointComputed || !storage_) {
    // Solve and store adjoint equation.
    for (unsigned i=0u; i<nobs_; ++i)
      M1_->applyInverse(*(p.get(i)),*(Mu.get(i)),z);
    p.scale(static_cast<Real>(-2));
    if (storage_) adjointStore_->set(p,param);
  }
}

template<typename Real>
ObjectivePV<Real>::ObjectivePV(const Ptr<MomentOperator<Real>>& M0,
                               const Ptr<MomentOperator<Real>>& M1,
                               const Ptr<Factors<Real>>& F,
                               bool storage)
  : M0_(M0), M1_(M1), F_(F),
    ts_(makePtr<TraceSampler<Real>>(F->createObservationVector(false))),
    stateStore_(makePtr<VectorController<Real>>()),
    MstateStore_(makePtr<VectorController<Real>>()),
    adjointStore_(makePtr<VectorController<Real>>()),
    rhsStore_(makePtr<VectorController<Real>>()),
    s_(F->createParameterVector(false)), q_(s_->clone()),
    e_(F->createObservationVector(false)),
    r1_(F->createParameterVector(true)), r2_(r1_->clone()),
    storage_(storage), nobs_(F->numObservations()) {
  std::vector<Ptr<Vector<Real>>> uvecs(nobs_);
  std::vector<Ptr<Vector<Real>>> Muvecs(nobs_);
  std::vector<Ptr<Vector<Real>>> pvecs(nobs_);
  std::vector<Ptr<Vector<Real>>> rvecs(nobs_);
  for (unsigned i=0; i<nobs_; ++i) {
    uvecs[i] = F->createParameterVector(false);
    Muvecs[i] = F->createParameterVector(true);
    pvecs[i] = F->createParameterVector(false);
    rvecs[i] = F->createParameterVector(true);
  }
  u_ = makePtr<PartitionedVector<Real>>(uvecs);
  Mu_ = makePtr<PartitionedVector<Real>>(Muvecs);
  p_ = makePtr<PartitionedVector<Real>>(pvecs);
  r_ = makePtr<PartitionedVector<Real>>(rvecs);
}

template<typename Real>
void ObjectivePV<Real>::update( const Vector<Real>& z, UpdateType type, int iter) {
  M0_->update(z,type,iter);
  M1_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
  MstateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectivePV<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  solveState(*u_,*Mu_,z,Objective<Real>::getParameter());
  // Assemble value
  Real val(0);
  for (unsigned i=0u; i<nobs_; ++i)
    val += (Mu_->get(i))->apply(*(u_->get(i)));
  return val;
}

template<typename Real>
void ObjectivePV<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();  
  // Solve state equation
  solveState(*u_,*Mu_,z,Objective<Real>::getParameter());
  // Solve adjoint equation
  solveAdjoint(*p_,*Mu_,z,Objective<Real>::getParameter());
  // Adjoint gradient
  g.zero();
  for (unsigned i=0u; i<nobs_; ++i) {
    M0_->applySampleMatrices(*g_,*(u_->get(i)),*(u_->get(i)));
    g.plus(*g_);
    M1_->applySampleMatrices(*g_,*(u_->get(i)),*(p_->get(i)));
    g.plus(*g_);
  }
}

template<typename Real>
void ObjectivePV<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  const Real one(1), two(2);
  if (g_ == nullPtr) g_ = hv.clone();  
  // Solve state equation
  solveState(*u_,*Mu_,z,Objective<Real>::getParameter());
  // Solve adjoint equation
  solveAdjoint(*p_,*Mu_,z,Objective<Real>::getParameter());
  hv.zero();
  for (unsigned i=0u; i<nobs_; ++i) {
    // Solve state sensitivity equation
    M1_->applyDeriv(*r1_,*(u_->get(i)),v);
    M1_->applyInverse(*s_,*r1_,z);
    // Solve adjoint sensitivity equation
    M0_->applyDeriv(*r1_,*(u_->get(i)),v);
    M1_->applyDeriv(*r2_,*(p_->get(i)),v);
    r2_->axpy(two,*r1_);
    M0_->apply(*r1_,*s_,z);
    r1_->scale(two);
    r1_->axpy(-one,*r2_);
    M1_->applyInverse(*q_,*r1_,z);
    // Assemble hessVec
    M1_->applySampleMatrices(*g_,*(u_->get(i)),*q_);
    hv.plus(*g_);
    M0_->applySampleMatrices(*g_,*(u_->get(i)),*s_);
    hv.axpy(-two,*g_);
    M1_->applySampleMatrices(*g_,*(p_->get(i)),*s_);
    hv.axpy(-one,*g_);
  }
}

} // End ROL::OED::Het Namespace

#endif
