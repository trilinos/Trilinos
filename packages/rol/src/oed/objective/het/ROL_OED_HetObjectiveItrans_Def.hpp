// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEITRANS_DEF_HPP
#define ROL_OED_HETOBJECTIVEITRANS_DEF_HPP

namespace ROL::OED::Het {

template<typename Real>
void ObjectiveItrans<Real>::solveState(Vector<Real>& u0, Vector<Real>& u1, Vector<Real>& Mu0, Vector<Real>& Mu1, const Vector<Real> &z, unsigned i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? state0Store_->get(u0,i) : false;
  state1Store_->get(u1,i);
  Mstate0Store_->get(Mu0,i);
  Mstate1Store_->get(Mu1,i);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    ts_->get(*r_,{static_cast<Real>(i)});
    M1_->applyInverse(u0,*r_,z);
    M1_->applyInverse(u1,*b_[i],z);
    M0_->apply(Mu0,u0,z);
    M0_->apply(Mu1,u1,z);
    // Store state.
    if (storage_) {
      state0Store_->set(u0,i);
      state1Store_->set(u1,i);
      Mstate0Store_->set(Mu0,i);
      Mstate1Store_->set(Mu1,i);
    }
  }
}

template<typename Real>
void ObjectiveItrans<Real>::solveAdjoint(Vector<Real>& p0, Vector<Real>& p1, const Vector<Real>& Mu0, const Vector<Real>& Mu1, const Vector<Real> &z, unsigned i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? adjoint0Store_->get(p0,i) : false;
  adjoint1Store_->get(p1,i);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M1_->applyInverse(p0,Mu1,z);
    M1_->applyInverse(p1,Mu0,z);
    p0.scale(static_cast<Real>(-1));
    p1.scale(static_cast<Real>(-1));
    // Store state.
    if (storage_) {
      adjoint0Store_->set(p0,i);
      adjoint1Store_->set(p1,i);
    }
  }
}

template<typename Real>
ObjectiveItrans<Real>::ObjectiveItrans( const Ptr<MomentOperator<Real>>& M0,
                                        const Ptr<MomentOperator<Real>>& M1,
                                        const Ptr<Factors<Real>>& F,
                                        const Ptr<SampleGenerator<Real>>& sampler,
                                        const Ptr<TraceSampler<Real>>& ts,
                                        const std::vector<Real>& wt,
                                        bool storage)
  : M0_(M0), M1_(M1), ts_(ts),
    state0Store_(makePtr<VectorController<Real,unsigned>>()),
    state1Store_(makePtr<VectorController<Real,unsigned>>()),
    Mstate0Store_(makePtr<VectorController<Real,unsigned>>()),
    Mstate1Store_(makePtr<VectorController<Real,unsigned>>()),
    adjoint0Store_(makePtr<VectorController<Real,unsigned>>()),
    adjoint1Store_(makePtr<VectorController<Real,unsigned>>()),
    u0_(F->createParameterVector(false)), u1_(u0_->clone()),
    Mu0_(u0_->dual().clone()), Mu1_(Mu0_->clone()),
    p0_(u0_->clone()), p1_(u1_->clone()), s_(u0_->clone()),
    r_(Mu0_->clone()), r1_(r_->clone()),
    storage_(storage), dim_(wt.size()) {
  auto obs = F->createObservationVector(false);
  b_.clear(); b_.resize(dim_);
  for (unsigned i = 0; i < dim_; ++i) {
    r_->zero();
    b_[i] = r_->clone();
    ts->get(*Mu0_,{static_cast<Real>(i)});
    for (int j = 0; j < sampler->numMySamples(); ++j) {
      F->apply(*obs,Mu0_->dual(),sampler->getMyPoint(j));
      F->applyAdjoint(*r1_,obs->dual(),sampler->getMyPoint(j));
      r_->axpy(sampler->getMyWeight(j),*r1_);
    }
    sampler->sumAll(*r_,*b_[i]);
    b_[i]->scale(wt[i]);
  }
}

template<typename Real>
void ObjectiveItrans<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  M0_->update(z,type,iter);
  M1_->update(z,type,iter);
  state0Store_->objectiveUpdate(type);
  state1Store_->objectiveUpdate(type);
  Mstate0Store_->objectiveUpdate(type);
  Mstate1Store_->objectiveUpdate(type);
  adjoint0Store_->objectiveUpdate(type);
  adjoint1Store_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectiveItrans<Real>::value( const Vector<Real> &z, Real &tol ) {
  Real val(0);
  for (unsigned i=0u; i<dim_; ++i) {
    // Solve state equation
    solveState(*u0_,*u1_,*Mu0_,*Mu1_,z,i);
    // Assemble value
    val += Mu1_->apply(*u0_);
  }
  return val;
}

template<typename Real>
void ObjectiveItrans<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = g.clone();
  g.zero();
  for (unsigned i=0u; i<dim_; ++i) {
    // Solve state equation
    solveState(*u0_,*u1_,*Mu0_,*Mu1_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p0_,*p1_,*Mu0_,*Mu1_,z,i);
    // Build gradient
    M0_->applySampleMatrices(*g_,*u0_,*u1_);
    g.plus(*g_);
    M1_->applySampleMatrices(*g_,*u0_,*p0_);
    g.plus(*g_);
    M1_->applySampleMatrices(*g_,*u1_,*p1_);
    g.plus(*g_);
  }
}

template<typename Real>
void ObjectiveItrans<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_==nullPtr) g_ = hv.clone();
  hv.zero();
  for (unsigned i = 0; i < dim_; ++i) {
    // Solve state equation
    solveState(*u0_,*u1_,*Mu0_,*Mu1_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p0_,*p1_,*Mu0_,*Mu1_,z,i);
    // Solve state sensitivity equation
    M1_->applyDeriv(*r_,*u0_,v);
    M1_->applyInverse(*s_,*r_,z);
    M1_->applySampleMatrices(*g_,*s_,*p0_);
    hv.axpy(static_cast<Real>(-1),*g_);
    M0_->applySampleMatrices(*g_,*s_,*u1_);
    hv.axpy(static_cast<Real>(-1),*g_);
    // Solve adjoint sensitivity equation
    M0_->apply(*r_,*s_,z);
    M0_->applyDeriv(*r1_,*u0_,v);
    r_->axpy(static_cast<Real>(-1),*r1_);
    M1_->applyDeriv(*r1_,*p1_,v);
    r_->axpy(static_cast<Real>(-1),*r1_);
    M1_->applyInverse(*s_,*r_,z);
    M1_->applySampleMatrices(*g_,*s_,*u1_);
    hv.plus(*g_);
    // Solve state sensitivity equation
    M1_->applyDeriv(*r_,*u1_,v);
    M1_->applyInverse(*s_,*r_,z);
    M1_->applySampleMatrices(*g_,*s_,*p1_);
    hv.axpy(static_cast<Real>(-1),*g_);
    M0_->applySampleMatrices(*g_,*s_,*u0_);
    hv.axpy(static_cast<Real>(-1),*g_);
    // Solve adjoint sensitivity equation
    M0_->apply(*r_,*s_,z);
    M0_->applyDeriv(*r1_,*u1_,v);
    r_->axpy(static_cast<Real>(-1),*r1_);
    M1_->applyDeriv(*r1_,*p0_,v);
    r_->axpy(static_cast<Real>(-1),*r1_);
    M1_->applyInverse(*s_,*r_,z);
    M1_->applySampleMatrices(*g_,*s_,*u0_);
    hv.plus(*g_);
  }
}

} // END ROL::OED::Het Namespace

#endif
