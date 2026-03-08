// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SEPARABLE_CONSTRAINT_DEF_H
#define ROL_SEPARABLE_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
SeparableConstraint<Real>::SeparableConstraint(const std::vector<Ptr<Constraint<Real>>> &cvec)
  : con_(cvec), size_(cvec.size()) {}

template<typename Real>
Ptr<Constraint<Real>> SeparableConstraint<Real>::get(unsigned ind) const {
  ROL_TEST_FOR_EXCEPTION( (size_ <= ind ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, get): "
                              "Dimension mismatch."); 
  return con_[ind];
}

template<typename Real>
void SeparableConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, update): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->update(*xp.get(i),type,iter);
}

template<typename Real>
void SeparableConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, update): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->update(*xp.get(i),flag,iter);
}

template<typename Real>
void SeparableConstraint<Real>::value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real>       &cp = dynamic_cast<PartitionedVector<Real>&>(c);
  const PartitionedVector<Real> &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != cp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, value): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, value): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->value(*cp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableConstraint<Real>::applyJacobian( Vector<Real> &jv,
                                         const Vector<Real> &v,
                                         const Vector<Real> &x,
                                         Real &tol ) {
  PartitionedVector<Real>       &jvp = dynamic_cast<PartitionedVector<Real>&>(jv);
  const PartitionedVector<Real>  &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != jvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyJacobian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyJacobian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyJacobian): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->applyJacobian(*jvp.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableConstraint<Real>::applyAdjointJacobian( Vector<Real> &ajv,
                                                const Vector<Real> &v,
                                                const Vector<Real> &x,
                                                Real &tol ) {
  PartitionedVector<Real>       &ajvp = dynamic_cast<PartitionedVector<Real>&>(ajv);
  const PartitionedVector<Real>   &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>   &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != ajvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointJacobian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointJacobian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointJacobian): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->applyAdjointJacobian(*ajvp.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableConstraint<Real>::applyAdjointHessian( Vector<Real> &ahuv,
                                               const Vector<Real> &u,
                                               const Vector<Real> &v,
                                               const Vector<Real> &x,
                                               Real &tol ) {
  PartitionedVector<Real>       &ahuvp = dynamic_cast<PartitionedVector<Real>&>(ahuv);
  const PartitionedVector<Real>    &up = dynamic_cast<const PartitionedVector<Real>&>(u);
  const PartitionedVector<Real>    &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>    &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != ahuvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointHessian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != up.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointHessian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointHessian): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyAdjointHessian): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->applyAdjointHessian(*ahuvp.get(i),*up.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableConstraint<Real>::applyPreconditioner(Vector<Real> &pv,
                                              const Vector<Real> &v,
                                              const Vector<Real> &x,
                                              const Vector<Real> &g,
                                              Real &tol) {
  PartitionedVector<Real>       &pvp = dynamic_cast<PartitionedVector<Real>&>(pv);
  const PartitionedVector<Real>  &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  const PartitionedVector<Real>  &gp = dynamic_cast<const PartitionedVector<Real>&>(g);
  ROL_TEST_FOR_EXCEPTION( (size_ != pvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyPreconditioner): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyPreconditioner): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyPreconditioner): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != gp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableConstraint, applyPreconditioner): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) con_[i]->applyPreconditioner(*pvp.get(i),*vp.get(i),*xp.get(i),*gp.get(i),tol);
}

template<typename Real>
void SeparableConstraint<Real>::setParameter(const std::vector<Real> &param) {
  Constraint<Real>::setParameter(param);
  for(unsigned i = 0; i < size_; ++i) con_[i]->setParameter(param);
}

} // namespace ROL

#endif
