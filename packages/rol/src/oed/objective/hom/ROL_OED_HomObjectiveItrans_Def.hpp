// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEITRANS_DEF_HPP
#define ROL_OED_HOMOBJECTIVEITRANS_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
void ObjectiveItrans<Real>::solveState(Vector<Real>& u, const Vector<Real> &z, int i) {
  std::vector<Real> param = {static_cast<Real>(i)};
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(u,i) : false;
  ts_->get(*r_,param);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M_->applyInverse(u,*r_,z);
    // Store state.
    if (storage_) stateStore_->set(u,i);
  }
}

template<typename Real>
void ObjectiveItrans<Real>::solveAdjoint(Vector<Real>& p, const Vector<Real> &z, int i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? adjointStore_->get(p,i) : false;
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M_->applyInverse(p,*b_[i],z);
    // Store state.
    if (storage_) adjointStore_->set(p,i);
  }
}

template<typename Real>
ObjectiveItrans<Real>::ObjectiveItrans( const Ptr<MomentOperator<Real>>& M,
                                        const Ptr<Factors<Real>>& F,
                                        const Ptr<SampleGenerator<Real>>& sampler,
                                        const Ptr<TraceSampler<Real>>& ts,
                                        const std::vector<Real>& wt,
                                        bool storage)
  : M_(M), ts_(ts),
    stateStore_(makePtr<VectorController<Real,int>>()),
    adjointStore_(makePtr<VectorController<Real,int>>()),
    u_(F->createParameterVector(false)), p_(u_->clone()), s_(u_->clone()),
    r_(u_->dual().clone()), storage_(storage), dim_(wt.size()) {
  int dim = wt.size();
  auto g0  = F->createParameterVector(true);
  auto b   = F->createParameterVector(true);
  auto obs = F->createObservationVector(false);
  b_.clear(); b_.resize(dim);
  for (int i = 0; i < dim; ++i) {
    b->zero();
    b_[i] = b->clone();
    ts->get(*g0,{static_cast<Real>(i)});
    for (int j = 0; j < sampler->numMySamples(); ++j) {
      F->apply(*obs,g0->dual(),sampler->getMyPoint(j));
      F->applyAdjoint(*r_,obs->dual(),sampler->getMyPoint(j));
      b->axpy(sampler->getMyWeight(j),*r_);
    }
    sampler->sumAll(*b,*b_[i]);
    b_[i]->scale(wt[i]);
  }
}

template<typename Real>
void ObjectiveItrans<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  M_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectiveItrans<Real>::value( const Vector<Real> &z, Real &tol ) {
  Real val(0);
  for (int i=0; i<dim_; ++i) {
    // Solve state equation
    solveState(*u_,z,i);
    // Assemble value
    val += b_[i]->apply(*u_);
  }
  return val;
}

template<typename Real>
void ObjectiveItrans<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = g.clone();
  g.zero();
  for (int i = 0; i < dim_; ++i) {
    // Solve state equation
    solveState(*u_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p_,z,i);
    // Build gradient
    M_->applySampleMatrices(*g_,*u_,*p_);
    g.plus(*g_);
  }
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void ObjectiveItrans<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = hv.clone();
  hv.zero();
  for (int i = 0; i < dim_; ++i) {
    // Solve state equation
    solveState(*u_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p_,z,i);
    // Solve state sensitivity equation
    M_->applyDeriv(*r_,*u_,v);
    M_->applyInverse(*s_,*r_,z);
    M_->applySampleMatrices(*g_,*s_,*p_);
    hv.plus(*g_);
    // Solve adjoint sensitivity equation
    M_->applyDeriv(*r_,*p_,v);
    M_->applyInverse(*s_,*r_,z);
    M_->applySampleMatrices(*g_,*s_,*u_);
    hv.plus(*g_);
  }
}

} // END ROL::OED::Hom Namespace

#endif
