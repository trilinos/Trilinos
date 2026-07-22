// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEA_DEF_HPP
#define ROL_OED_HETOBJECTIVEA_DEF_HPP

namespace ROL::OED::Het {

template<typename Real>
void ObjectiveA<Real>::solveState(Vector<Real>& u, Vector<Real>& Mu, const Vector<Real>& z, unsigned i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(u,i) : false;
  MstateStore_->get(Mu,i);
  if (!isComputed || !storage_) {
    std::vector<Real> param = {static_cast<Real>(i)};
    ts_->get(*r1_,param);
    M1_->applyInverse(u,*r1_,z);
    M0_->apply(Mu,u,z);
    if (storage_) {
      stateStore_->set(u,i);
      MstateStore_->set(Mu,i);
    }
  }
}

template<typename Real>
void ObjectiveA<Real>::solveAdjoint(Vector<Real>& p, const Vector<Real>& Mu, const Vector<Real>& z, unsigned i) {
  // Check if state has been computed.
  bool isComputed = storage_ ? adjointStore_->get(p,i) : false;
  if (!isComputed || !storage_) {
    M1_->applyInverse(p,Mu,z);
    p.scale(static_cast<Real>(-2));
    if (storage_) adjointStore_->set(p,i);
  }
}

template<typename Real>
ObjectiveA<Real>::ObjectiveA(const Ptr<MomentOperator<Real>>& M0,
                             const Ptr<MomentOperator<Real>>& M1,
                             const Ptr<const Vector<Real>>& theta,
                             const Ptr<TraceSampler<Real>>& ts,
                             const std::vector<Real>& weight,
                             bool storage)
  : M0_(M0), M1_(M1), ts_(ts), weight_(weight),
    stateStore_(makePtr<VectorController<Real,unsigned>>()),
    MstateStore_(makePtr<VectorController<Real,unsigned>>()),
    adjointStore_(makePtr<VectorController<Real,unsigned>>()),
    u_(theta->clone()), Mu_(theta->dual().clone()), p_(u_->clone()),
    r1_(Mu_->clone()), r2_(Mu_->clone()),
    s_(u_->clone()), q_(p_->clone()),
    storage_(storage), dim_(weight.size()) {}

template<typename Real>
void ObjectiveA<Real>::update( const Vector<Real>& z, UpdateType type, int iter) {
  M0_->update(z,type,iter);
  M1_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
  MstateStore_->objectiveUpdate(type);
  adjointStore_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectiveA<Real>::value( const Vector<Real> &z, Real &tol ) {
  Real val(0);
  for (unsigned i=0; i<dim_; ++i) {
    // Solve state equation
    solveState(*u_,*Mu_,z,i);
    // Assemble objective function value
    val += weight_[i] * Mu_->apply(*u_);
  }
  return val;
}

template<typename Real>
void ObjectiveA<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();  
  g.zero();
  for (unsigned i=0; i<dim_; ++i) {
    // Solve state equation
    solveState(*u_,*Mu_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p_,*Mu_,z,i);
    // Adjoint gradient
    M0_->applySampleMatrices(*g_,*u_,*u_);
    g.axpy(weight_[i],*g_);
    M1_->applySampleMatrices(*g_,*u_,*p_);
    g.axpy(weight_[i],*g_);
  }
}

template<typename Real>
void ObjectiveA<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = hv.clone();  
  hv.zero();
  const Real one(1), two(2);
  for (unsigned i=0; i<dim_; ++i) {
    // Solve state equation
    solveState(*u_,*Mu_,z,i);
    // Solve adjoint equation
    solveAdjoint(*p_,*Mu_,z,i);
    // Solve state sensitivity equation
    M1_->applyDeriv(*r1_,*u_,v);
    M1_->applyInverse(*s_,*r1_,z);
    // Solve adjoint sensitivity equation
    M0_->applyDeriv(*r1_,*u_,v);
    M1_->applyDeriv(*r2_,*p_,v);
    r2_->axpy(two,*r1_);
    M0_->apply(*r1_,*s_,z);
    r1_->scale(two);
    r1_->axpy(-one,*r2_);
    M1_->applyInverse(*q_,*r1_,z);
    // Assemble hessVec
    M1_->applySampleMatrices(*g_,*u_,*q_);
    hv.axpy(weight_[i],*g_);
    M0_->applySampleMatrices(*g_,*u_,*s_);
    hv.axpy(-two*weight_[i],*g_);
    M1_->applySampleMatrices(*g_,*p_,*s_);
    hv.axpy(-one*weight_[i],*g_);
  }
}

} // End ROL::OED::Het Namespace

#endif
