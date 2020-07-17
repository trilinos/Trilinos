// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
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
LocalCrsMatrixOperator (const std::shared_ptr<local_matrix_type>& A)
  : A_ (A)
{
  const char tfecfFuncName[] = "LocalCrsMatrixOperator: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (A_.get () == nullptr, std::invalid_argument,
     "Input matrix A is null.");
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
  //Only create A_ordinal_rowptrs if:
  //  - KokkosKernels cuSPARSE support is enabled (otherwise, no benefit)
  //  - The execution space is CUDA
  //  - The local matrix offset and ordinal types are different (otherwise, no reason to enable)
  //  - The number of entries can be represented by the ordinal type.
  using kk_offset_t = typename std::remove_const<typename local_matrix_type::size_type>::type;
  using kk_ordinal_t = typename std::remove_const<typename local_matrix_type::ordinal_type>::type;
  using exec_space = typename Device::execution_space;
  if(std::is_same<exec_space, Kokkos::Cuda>::value &&
      !std::is_same<kk_offset_t, kk_ordinal_t>::value &&
      A_->nnz() < static_cast<kk_offset_t>(Teuchos::OrdinalTraits<kk_ordinal_t>::max()))
  {
    A_ordinal_rowptrs = ordinal_view_type(Kokkos::ViewAllocateWithoutInitializing("A_ordinal_rowptrs"), A_->numRows() + 1);
    //This is just like a deep copy, but it implicitly converts each element
    KokkosKernels::Impl::copy_view<typename local_graph_type::row_map_type, ordinal_view_type, exec_space>
      (A_ordinal_rowptrs.extent(0), A_->graph.row_map, A_ordinal_rowptrs);
    A_cusparse = local_cusparse_matrix_type("A(cusparse)", A_->numRows(), A_->numCols(), A_->nnz(), A_->values, A_ordinal_rowptrs, A_->graph.entries);
  }
#endif
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
  //Currently KK has no cusparse wrapper for rank-2 (SpMM)
  //TODO: whent that is supported, use A_cusparse for that case also
  if(X.extent(1) == size_t(1) && A_ordinal_rowptrs.extent(0))
  {
    KokkosSparse::spmv (op, alpha, A_cusparse, Kokkos::subview(X, Kokkos::ALL(), 0),
                            beta, Kokkos::subview(Y, Kokkos::ALL(), 0));
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
  const bool conjugate = (mode == Teuchos::CONJ_TRANS);
  const bool transpose = (mode != Teuchos::NO_TRANS);

#ifdef HAVE_TPETRA_DEBUG
  const char tfecfFuncName[] = "applyLoadBalanced: ";

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
  //Select the merge path algorithm (used if available, otherwise has no effect)
  //TODO BMK: If/when KokkosKernels gets its own SPMV implementation for imbalanced rows,
  //call that here or select it using Controls.
  //Ideally it supports multivectors from the beginning.
  if((Details::Behavior::useMergePathMultiVector() || X.extent(1) == size_t(1)) && A_ordinal_rowptrs.extent(0))
  {
    KokkosKernels::Experimental::Controls controls;
    controls.setParameter("algorithm", "merge");
    //Apply on one column at a time (must be rank-1)
    for(size_t vec = 0; vec < X.extent(1); vec++)
    {
      KokkosSparse::spmv (controls, op,
          alpha, A_cusparse, Kokkos::subview(X, Kokkos::ALL(), vec),
          beta, Kokkos::subview(Y, Kokkos::ALL(), vec));
    }
  }
  else
  {
    //Just run multivector version of spmv (no controls, and no cusparse support)
    KokkosSparse::spmv (op, alpha, *A_, X, beta, Y);
  }
}

template<class MultiVectorScalar, class MatrixScalar, class Device>
const typename LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::local_matrix_type&
LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>::
getLocalMatrix () const
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

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_DEF_HPP
