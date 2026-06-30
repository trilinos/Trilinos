// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLOCKDIAGONALOPERATOR_DEF_H
#define ROL_BLOCKDIAGONALOPERATOR_DEF_H

#include "ROL_PartitionedVector.hpp"

namespace ROL {

template<class Real>
BlockDiagonalOperator<Real>::BlockDiagonalOperator( const std::vector<Ptr<LinearOperator<Real>>> &diag )
  : diag_(diag), numOp_(diag.size()) {}

template<class Real>
const Ptr<LinearOperator<Real>> BlockDiagonalOperator<Real>::get(unsigned i) const {
  ROL_TEST_FOR_EXCEPTION( (i >= numOp_), std::invalid_argument,
                              ">>> ERROR (ROL_BlockDiagonalOperator, get): "
                              "Input exceeds number of blocks.");
  return diag_[i];
}

template<class Real>
void BlockDiagonalOperator<Real>::apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
  // Downcast to Partitioned Vectors
  PartitionedVector<Real>       &Hv_part = dynamic_cast<PartitionedVector<Real>&>(Hv);
  const PartitionedVector<Real>  &v_part = dynamic_cast<const PartitionedVector<Real>&>(v);

  unsigned nvec1 = v_part.numVectors();
  unsigned nvec2 = Hv_part.numVectors();

  ROL_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                              ">>> ERROR (ROL_BlockDiagonalOperator, apply): "
                              "Mismatch between input and output number of subvectors.");

  ROL_TEST_FOR_EXCEPTION( (numOp_ != nvec2 ) , std::invalid_argument, 
                              ">>> ERROR (ROL_BlockDiagonalOperator, apply): "
                              "Block operator dimension mismatch."); 

  for (unsigned i=0; i<numOp_; ++i)
    diag_[i]->apply(*(Hv_part.get(i)), *(v_part.get(i)), tol);
}

template<class Real>
void BlockDiagonalOperator<Real>::applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
  // Downcast to Partitioned Vectors
  PartitionedVector<Real>       &Hv_part = dynamic_cast<PartitionedVector<Real>&>(Hv);
  const PartitionedVector<Real>  &v_part = dynamic_cast<const PartitionedVector<Real>&>(v);

  unsigned nvec1 = v_part.numVectors();
  unsigned nvec2 = Hv_part.numVectors();

  ROL_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyInverse): "
                              "Mismatch between input and output number of subvectors.");

  ROL_TEST_FOR_EXCEPTION( (numOp_ != nvec2 ) , std::invalid_argument, 
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyInverse): "
                              "Block operator dimension mismatch."); 

  for (unsigned i=0; i<numOp_; ++i)
    diag_[i]->applyInverse(*(Hv_part.get(i)), *(v_part.get(i)), tol);
}

template<class Real>
void BlockDiagonalOperator<Real>::applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
  // Downcast to Partitioned Vectors
  PartitionedVector<Real>       &Hv_part = dynamic_cast<PartitionedVector<Real>&>(Hv);
  const PartitionedVector<Real>  &v_part = dynamic_cast<const PartitionedVector<Real>&>(v);

  unsigned nvec1 = v_part.numVectors();
  unsigned nvec2 = Hv_part.numVectors();

  ROL_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyAdjoint): "
                              "Mismatch between input and output number of subvectors.");

  ROL_TEST_FOR_EXCEPTION( (numOp_ != nvec2 ) , std::invalid_argument, 
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyAdjoint): "
                              "Block operator dimension mismatch."); 

  for (unsigned i=0; i<numOp_; ++i)
    diag_[i]->applyAdjoint(*(Hv_part.get(i)), *(v_part.get(i)), tol);
}

template<class Real>
void BlockDiagonalOperator<Real>::applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
  // Downcast to Partitioned Vectors
  PartitionedVector<Real>       &Hv_part = dynamic_cast<PartitionedVector<Real>&>(Hv);
  const PartitionedVector<Real>  &v_part = dynamic_cast<const PartitionedVector<Real>&>(v);

  unsigned nvec1 = v_part.numVectors();
  unsigned nvec2 = Hv_part.numVectors();

  ROL_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyAdjointInverse): "
                              "Mismatch between input and output number of subvectors.");

  ROL_TEST_FOR_EXCEPTION( (numOp_ != nvec2 ) , std::invalid_argument, 
                              ">>> ERROR (ROL_BlockDiagonalOperator, applyAdjointInverse): "
                              "Block operator dimension mismatch."); 

  for (unsigned i=0; i<numOp_; ++i)
    diag_[i]->applyAdjointInverse(*(Hv_part.get(i)), *(v_part.get(i)), tol);
}

} // namespace ROL

#endif // ROL_BLOCKDIAGONALOPERATOR_DEF_H
