// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_COMPOSITECONSTRAINT_SIMOPT_DEF_H
#define ROL_COMPOSITECONSTRAINT_SIMOPT_DEF_H

namespace ROL {

template<typename Real>
CompositeConstraint_SimOpt<Real>::CompositeConstraint_SimOpt(
                                  const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
                                  const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
                                  const Vector<Real> &cVal,
                                  const Vector<Real> &cRed,
                                  const Vector<Real> &u,
                                  const Vector<Real> &Sz,
                                  const Vector<Real> &z,
                                  bool storage,
                                  bool isConRedParametrized)
  : Constraint_SimOpt<Real>(), conVal_(conVal), conRed_(conRed),
    updateFlag_(true), newUpdate_(false), updateIter_(0), updateType_(UpdateType::Initial),
    storage_(storage), isConRedParametrized_(isConRedParametrized) {
  Sz_      = Sz.clone();
  primRed_ = cRed.clone();
  dualRed_ = cRed.dual().clone();
  primZ_   = z.clone();
  dualZ_   = z.dual().clone();
  dualZ1_  = z.dual().clone();
  primU_   = u.clone();
  stateStore_ = ROL::makePtr<VectorController<Real>>();
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update(const Vector<Real> &u,
                                              const Vector<Real> &z,
                                              bool flag, int iter) {
  update_1(u,flag,iter);
  update_2(z,flag,iter);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update_1(const Vector<Real> &u,
                                                bool flag, int iter) {
  primU_->set(u);
  conVal_->update_1(u, flag, iter);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update_2(const Vector<Real> &z,
                                                bool flag, int iter) {
  conRed_->update_2(z,flag,iter);
  updateFlag_ = flag;
  updateIter_ = iter;
  stateStore_->constraintUpdate(true);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update(const Vector<Real> &u,
                                              const Vector<Real> &z,
                                              UpdateType type, int iter) {
  update_1(u,type,iter);
  update_2(z,type,iter);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update_1(const Vector<Real> &u,
                                                UpdateType type, int iter) {
  primU_->set(u);
  conVal_->update_1(u,type,iter);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::update_2(const Vector<Real> &z,
                                                UpdateType type, int iter) {
  conRed_->update_2(z,type,iter);
  updateType_ = type;
  updateIter_ = iter;
  stateStore_->constraintUpdate(type);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::solve_update(const Vector<Real> &u,
                                                    const Vector<Real> &z,
                                                    UpdateType type, int iter) {
  Real ctol(std::sqrt(ROL_EPSILON<Real>()));
  solveConRed(*Sz_, z, ctol);
  conVal_->solve_update(u,*Sz_,type,iter);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::value(Vector<Real> &c,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->value(c, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::solve(Vector<Real> &c,
                                             Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->solve(c, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyJacobian_1(Vector<Real> &jv,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyJacobian_1(jv, v, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyJacobian_2(Vector<Real> &jv,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &u,
                                                       const Vector<Real> &z,
                                                       Real &tol) { 
  solveConRed(*Sz_, z, tol);
  applySens(*primZ_, v, *Sz_, z, tol);
  conVal_->applyJacobian_2(jv, *primZ_, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyInverseJacobian_1(Vector<Real> &ijv,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyInverseJacobian_1(ijv, v, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointJacobian_1(Vector<Real> &ajv,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyAdjointJacobian_1(ajv, v, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointJacobian_2(Vector<Real> &ajv,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyAdjointJacobian_2(*dualZ_, v, u, *Sz_, tol);
  applyAdjointSens(ajv, *dualZ_, *Sz_, z, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyInverseAdjointJacobian_1(Vector<Real> &ijv,
                                                                     const Vector<Real> &v,
                                                                     const Vector<Real> &u,
                                                                     const Vector<Real> &z,
                                                                     Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyInverseAdjointJacobian_1(ijv, v, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointHessian_11(Vector<Real> &ahwv,
                                                              const Vector<Real> &w,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyAdjointHessian_11(ahwv, w, v, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointHessian_12(Vector<Real> &ahwv,
                                                              const Vector<Real> &w,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  conVal_->applyAdjointHessian_12(*dualZ_, w, v, u, *Sz_, tol);
  applyAdjointSens(ahwv, *dualZ_, *Sz_, z, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointHessian_21(Vector<Real> &ahwv,
                                                              const Vector<Real> &w,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  solveConRed(*Sz_, z, tol);
  applySens(*primZ_, v, *Sz_, z, tol);
  conVal_->applyAdjointHessian_21(ahwv, w, *primZ_, u, *Sz_, tol);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointHessian_22(Vector<Real> &ahwv,
                                                              const Vector<Real> &w,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &u,
                                                              const Vector<Real> &z,
                                                              Real &tol) {
  ahwv.zero();
  solveConRed(*Sz_, z, tol);
  applySens(*primZ_, v, *Sz_, z, tol);

  conVal_->applyAdjointJacobian_2(*dualZ_, w, u, *Sz_, tol);
  conRed_->applyInverseAdjointJacobian_1(*dualRed_, *dualZ_, *Sz_, z, tol);
  conRed_->applyAdjointHessian_22(*dualZ_, *dualRed_, v, *Sz_, z, tol);
  ahwv.axpy(static_cast<Real>(-1), *dualZ_);
  conRed_->applyAdjointHessian_12(*dualZ_, *dualRed_, *primZ_, *Sz_, z, tol);
  ahwv.axpy(static_cast<Real>(-1), *dualZ_);

  conRed_->applyAdjointHessian_11(*dualZ1_, *dualRed_, *primZ_, *Sz_, z, tol);
  conRed_->applyAdjointHessian_21(*dualZ_, *dualRed_, v, *Sz_, z, tol);
  dualZ1_->plus(*dualZ_); 
  dualZ1_->scale(static_cast<Real>(-1));
  
  conVal_->applyAdjointHessian_22(*dualZ_, w, *primZ_, u, *Sz_, tol);
  dualZ1_->plus(*dualZ_); 

  applyAdjointSens(*dualZ_, *dualZ1_, *Sz_, z, tol);
  ahwv.plus(*dualZ_);
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::setParameter(const std::vector<Real> &param) {
  conVal_->setParameter(param);
  if (isConRedParametrized_) {
    conRed_->setParameter(param);
    Constraint_SimOpt<Real>::setParameter(param);
  }
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::solveConRed(Vector<Real> &Sz,
                                                   const Vector<Real> &z,
                                                   Real &tol) {
  std::vector<Real> param = Constraint_SimOpt<Real>::getParameter();
  // Check if state has been computed.
  bool isComputed = false;
  if (storage_) isComputed = stateStore_->get(Sz,param);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    conRed_->solve(*primRed_,Sz,z,tol);
    // Update equality constraint with new Sim variable.
    if (newUpdate_) conRed_->update_1(Sz,updateType_,updateIter_);
    else            conRed_->update_1(Sz,updateFlag_,updateIter_);
    // Update equality constraint.
    if (newUpdate_) conRed_->update(Sz,z,updateType_,updateIter_);
    else            conRed_->update(Sz,z,updateFlag_,updateIter_);
    if (newUpdate_) conVal_->update(*primU_,Sz,updateType_,updateIter_);
     else           conVal_->update(*primU_,Sz,updateFlag_,updateIter_);
    // Store state.
    if (storage_) stateStore_->set(Sz,param);
  }
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applySens(Vector<Real> &jv,
                                                 const Vector<Real> &v,
                                                 const Vector<Real> &Sz,
                                                 const Vector<Real> &z,
                                                 Real &tol) { 
  // Solve linearization of reducible constraint in direction v
  conRed_->applyJacobian_2(*primRed_, v, Sz, z, tol);
  conRed_->applyInverseJacobian_1(jv, *primRed_, Sz, z, tol);
  jv.scale(static_cast<Real>(-1));
}

template<typename Real>
void CompositeConstraint_SimOpt<Real>::applyAdjointSens(Vector<Real> &ajv,
                                                        const Vector<Real> &v,
                                                        const Vector<Real> &Sz,
                                                        const Vector<Real> &z,
                                                        Real &tol) {
  // Solve adjoint of linearized reducible constraint
  conRed_->applyInverseAdjointJacobian_1(*dualRed_, v, Sz, z, tol);
  conRed_->applyAdjointJacobian_2(ajv, *dualRed_, Sz, z, tol);
  ajv.scale(static_cast<Real>(-1));
}

} // namespace ROL

#endif
