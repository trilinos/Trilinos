// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_LOCALCRSMATRIXOPERATOR_DEF_HPP
#define TPETRA_LOCALCRSMATRIXOPERATOR_DEF_HPP

#include "Tpetra_LocalOperator.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "KokkosSparse.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_OrdinalTraits.hpp"

namespace Tpetra {

template<class MultiVectorScalar, class MatrixScalar, class Device>
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_device_type>& A)
  : A_ (A), have_A_cusparse(false)
{
  const char tfecfFuncName[] = "LocalCrsMatrixOperator: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (A_.get () == nullptr, std::invalid_argument,
     "Input matrix A is null.");
}

template<class MultiVectorScalar, class MatrixScalar, class Device>
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_device_type>& A, const ordinal_view_type& A_ordinal_rowptrs) :
  A_ (A),
  A_cusparse("LocalCrsMatrixOperator_cuSPARSE", A->numRows(), A->numCols(), A->nnz(),
      A->values, A_ordinal_rowptrs, A->graph.entries),
  have_A_cusparse(true)
{
  const char tfecfFuncName[] = "LocalCrsMatrixOperator: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (A_.get () == nullptr, std::invalid_argument,
     "Input matrix A is null.");
}

template<class MultiVectorScalar, class MatrixScalar, class Device>
bool
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
hasTransposeApply () const
{
  return true;
}

template<class MultiVectorScalar, class MatrixScalar, class Device>
void
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
apply (Kokkos::View<const mv_scalar_type**, array_layout,
         device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
       Kokkos::View<mv_scalar_type**, array_layout,
         device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
       const Teuchos::ETransp mode,
       const mv_scalar_type alpha,
       const mv_scalar_type beta) const
{
  const bool conjugate = (mode == Teuchos::CONJ_TRANS);
  const bool transpose = (mode != Teuchos::NO_TRANS);

#ifdef HAVE_TPETRA_DEBUG
  const char tfecfFuncName[] = "apply: ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (X.extent (1) != Y.extent (1), std::runtime_error,
     "X.extent(1) = " << X.extent (1) << " != Y.extent(1) = "
     << Y.extent (1) << ".");
  // If the two pointers are NULL, then they don't alias one
  // another, even though they are equal.
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (X.data () == Y.data () && X.data () != nullptr,
     std::runtime_error, "X and Y may not alias one another.");
#endif // HAVE_TPETRA_DEBUG

  const auto op = transpose ?
    (conjugate ? KokkosSparse::ConjugateTranspose :
     KokkosSparse::Transpose) : KokkosSparse::NoTranspose;
  if(have_A_cusparse)
  {
    KokkosSparse::spmv (op, alpha, A_cusparse, X, beta, Y);
  }
  else
  {
    KokkosSparse::spmv (op, alpha, *A_, X, beta, Y);
  }
}

/// \brief Same behavior as \c apply() above, except give KokkosKernels a hint to use
///  an SPMV algorithm that can efficiently handle matrices with imbalanced rows.
template<class MultiVectorScalar, class MatrixScalar, class Device>
void
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
applyImbalancedRows (
       Kokkos::View<const mv_scalar_type**, array_layout,
         device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
       Kokkos::View<mv_scalar_type**, array_layout,
         device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
       const Teuchos::ETransp mode,
       const mv_scalar_type alpha,
       const mv_scalar_type beta) const
{
  apply(X, Y, mode, alpha, beta);
}

template<class MultiVectorScalar, class MatrixScalar, class Device>
const typename LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::local_matrix_device_type&
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
getLocalMatrixDevice () const
{
  return *A_;
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

// We only explicitly instantiate for MultiVectorScalar ==
// MatrixScalar, which is what CrsMatrix needs.

#define TPETRA_LOCALCRSMATRIXOPERATOR_INSTANT(SC,NT) \
  template class LocalCrsMatrixOperator< SC, SC, NT::device_type >;

// If we want mixed versions, we use this macro.

#define TPETRA_LOCALCRSMATRIXOPERATOR_MIXED_INSTANT(SC,MATSC,LO,GO,NT)    \
  template class LocalCrsMatrixOperator< SC, MATSC, NT::device_type >;

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_DEF_HPP
