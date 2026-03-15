// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SEPARABLEOBJECTIVE_DEF_H
#define ROL_SEPARABLEOBJECTIVE_DEF_H

namespace ROL {

template<typename Real>
SeparableObjective<Real>::SeparableObjective(const std::vector<Ptr<Objective<Real>>> &obj)
  : obj_(obj), size_(obj.size()) {}

template<typename Real>
const Ptr<Objective<Real>> SeparableObjective<Real>::get(unsigned i) const {
  ROL_TEST_FOR_EXCEPTION( (size_ <= i ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, get): "
                              "Dimension mismatch."); 
  return obj_[i];
}

template<typename Real>
void SeparableObjective<Real>::update(const Vector<Real> &x, UpdateType type, int iter) {
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, update): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->update(*xp.get(i),type,iter);
}

template<typename Real>
void SeparableObjective<Real>::update(const Vector<Real> &x, bool flag, int iter) {
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, update): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->update(*xp.get(i),flag,iter);
}

template<typename Real>
Real SeparableObjective<Real>::value( const Vector<Real> &x, Real &tol ) {
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, value): "
                              "Dimension mismatch."); 
  Real val(0);
  for(unsigned i = 0; i < size_; ++i) val += obj_[i]->value(*xp.get(i),tol);
  return val;
}

template<typename Real>
void SeparableObjective<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real>       &gp = dynamic_cast<PartitionedVector<Real>&>(g);
  const PartitionedVector<Real> &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != gp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, gradient): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, gradient): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->gradient(*gp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableObjective<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real>       &hvp = dynamic_cast<PartitionedVector<Real>&>(hv);
  const PartitionedVector<Real>  &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != hvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, gradient): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->hessVec(*hvp.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableObjective<Real>::invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real>       &hvp = dynamic_cast<PartitionedVector<Real>&>(hv);
  const PartitionedVector<Real>  &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != hvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, gradient): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->invHessVec(*hvp.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableObjective<Real>::precond( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  PartitionedVector<Real>       &hvp = dynamic_cast<PartitionedVector<Real>&>(hv);
  const PartitionedVector<Real>  &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
  const PartitionedVector<Real>  &xp = dynamic_cast<const PartitionedVector<Real>&>(x);
  ROL_TEST_FOR_EXCEPTION( (size_ != hvp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, gradient): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != vp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  ROL_TEST_FOR_EXCEPTION( (size_ != xp.numVectors() ) , std::invalid_argument, 
                              ">>> ERROR (ROL_SeparableObjective, hessVec): "
                              "Dimension mismatch."); 
  for(unsigned i = 0; i < size_; ++i) obj_[i]->precond(*hvp.get(i),*vp.get(i),*xp.get(i),tol);
}

template<typename Real>
void SeparableObjective<Real>::setParameter(const std::vector<Real> &param) {
  Objective<Real>::setParameter(param);
  for(unsigned i = 0; i < size_; ++i) obj_[i]->setParameter(param);
}

} // namespace ROL

#endif
