// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_CRSMATRIXASSEMBLEELEMENT_HPP
#define TPETRA_DETAILS_CRSMATRIXASSEMBLEELEMENT_HPP

#include "KokkosSparse_CrsMatrix.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {

/// \brief <tt> A(lclRow, lclColsInds[sortPerm[j]]) += vals[sortPerm[j]]</tt>,
///   for all j in </tt>0 .. eltDim-1</tt>.
///
/// In the row of the matrix A with the local row index lclRow, find
/// entries with column indices lclColInds, and sum into those entries
/// with vals.  Assume that lclColInds[sortPerm] is sorted, and that
/// the column indices in that row of the matrix are sorted as well.
/// Use linear search to find the entries in that row of the matrix.
///
/// \tparam SparseMatrixType Specialization of KokkosSparse::CrsMatrix.
/// \tparam ValsViewType Specialization of a 1-D Kokkos::View.
///
/// \param A [in/out] Sparse matrix whose entries to modify.
/// \param lclRow [in] Local index of the row in the matrix A to
///   modify.  lclRow MUST be a valid local row index of A.
/// \param lclColInds [in] Local column indices to modify in that row.
/// \param sortPerm [in] Permutation that makes lclColInds sorted.
///   That is, <tt>lclColInds[sortPerm]</tt> is sorted.
/// \param vals [in] Input 1-D Kokkos::View of the values to ruse.
///   This is a Kokkos::View and not a raw 1-D array, because it may
///   be strided, if the original element being used (see
///   crsMatrixSumInElement) has a column-major layout.
/// \param numEntInInput [in] Number of entries in the input.  This
///   function will read the first numEntInInput entries of
///   lclColInds, sortPerm, and vals.
/// \param forceAtomic [in] Whether to use atomic updates when
///   modifying the entries of the matrix A.  This MUST be a
///   compile-time constant.  It defaults to whether the matrix's
///   Kokkos execution space is NOT Kokkos::Serial.
/// \param checkInputIndices [in] Whether to check whether the input
///   indices are valid column indices before just using them.  For
///   forwards compatibility, this should always be a compile-time
///   constant.  Default is true, that is, always check.
///
/// \return If checkInputIndices is true, return the number of input
///   indices that are valid column indices in that row of the matrix.
///   If checkInputIndices is false, just return numEntInInput.
template<class SparseMatrixType,
         class ValsViewType>
KOKKOS_FUNCTION
typename SparseMatrixType::ordinal_type
crsMatrixSumIntoValues_sortedSortedLinear (const SparseMatrixType& A,
                                           const typename SparseMatrixType::ordinal_type lclRow,
                                           const typename SparseMatrixType::ordinal_type lclColInds[],
                                           const typename SparseMatrixType::ordinal_type sortPerm[],
                                           const ValsViewType& vals,
                                           const typename SparseMatrixType::ordinal_type numEntInInput,
                                           const bool forceAtomic =
#ifdef KOKKOS_ENABLE_SERIAL
                                           ! std::is_same<typename SparseMatrixType::device_type::execution_space, Kokkos::Serial>::type,
#else // NOT KOKKOS_ENABLE_SERIAL
                                           false,
#endif // KOKKOS_ENABLE_SERIAL
                                           const bool checkInputIndices = true)
{
  typedef typename std::remove_const<typename SparseMatrixType::value_type>::type
    matrix_scalar_type;
  static_assert (std::is_same<matrix_scalar_type,
                 typename SparseMatrixType::value_type>::value,
                 "The matrix's entries must have a nonconst type.");
  // static_assert (std::is_assignable<matrix_scalar_type,
  //                typename std::decay< decltype (A.values[0] + vals[0]) >::type>::value,
  //                "The result of adding a matrix entry and an entry of vals "
  //                "MUST be assignable to a matrix entry.");
  typedef typename SparseMatrixType::ordinal_type LO;
  static_assert (std::is_integral<LO>::value, "SparseMatrixType::ordinal_type "
                 "must be a built-in integer type.");

  // If lclRow is NOT a valid row index, this will return a view of
  // zero entries.  If checkInputIndices is true, thus, then none of
  // the input indices will be valid in that case.
  auto row_view = A.row (lclRow);
  const LO numEntInRow = static_cast<LO> (row_view.length);
  // Number of valid local column indices found, that is, the number
  // of input indices that are valid column indices found in row
  // lclRow of the matrix.  If not checking, we just return the number
  // of input indices.
  LO numValid = checkInputIndices ? static_cast<LO> (0) : numEntInRow;

  // Since both the matrix row and the input (after permutation) are
  // sorted, we only need to pass once over the matrix row.  'offset'
  // tells us the current search position in the matrix row.
  LO offset = 0;
  for (LO j = 0; j < numEntInInput; ++j) {
    const LO perm_index = sortPerm[j];
    const LO lclColInd = lclColInds[perm_index];
    // Search linearly in the matrix row for the current index.
    // If we ever want binary search, this would be the place.
    while (row_view.colidx(offset) != lclColInd) {
      ++offset;
    }

    // If we could make checkInputIndices a compile-time constant,
    // then the compiler might not need to insert a branch here.  This
    // should help vectorization, if vectorization is possible.
    if (checkInputIndices) {
      if (offset != numEntInRow) {
        // If we could make forceAtomic a compile-time constant, then
        // the compiler might not need to insert a branch here.  This
        // should help vectorization, if vectorization is possible.
        if (forceAtomic) {
          Kokkos::atomic_add (&(row_view.value(offset)), vals[perm_index]);
        }
        else {
          row_view.value(offset) += vals[perm_index];
        }
        ++numValid;
      }
    }
    else { // don't check input indices; assume they are in the row
      // See above note on forceAtomic.
      if (forceAtomic) {
        Kokkos::atomic_add (&(row_view.value(offset)), vals[perm_index]);
      }
      else {
        row_view.value(offset) += vals[perm_index];
      }
    }
  }

  return numValid;
}

/// \brief <tt> A(lclRow, lclColsInds[sortPerm[j]]) = vals[sortPerm[j]]</tt>,
///   for all j in </tt>0 .. eltDim-1</tt>.
///
/// In the row of the matrix A with the local row index lclRow, find
/// entries with column indices lclColInds, and replace those entries
/// with vals.  Assume that lclColInds[sortPerm] is sorted, and that
/// the column indices in that row of the matrix are sorted as well.
/// Use linear search to find the entries in that row of the matrix.
///
/// \tparam SparseMatrixType Specialization of KokkosSparse::CrsMatrix.
/// \tparam ValsViewType Specialization of a 1-D Kokkos::View.
///
/// \param A [in/out] Sparse matrix whose entries to modify.
/// \param lclRow [in] Local index of the row in the matrix A to
///   modify.  lclRow MUST be a valid local row index of A.
/// \param lclColInds [in] Local column indices to modify in that row.
/// \param sortPerm [in] Permutation that makes lclColInds sorted.
///   That is, <tt>lclColInds[sortPerm]</tt> is sorted.
/// \param vals [in] Input 1-D Kokkos::View of the values to use.
///   This is a Kokkos::View and not a raw 1-D array, because it may
///   be strided, if the original element being used (see
///   crsMatrixSumInElement) has a column-major layout.
/// \param numEntInInput [in] Number of entries in the input.  This
///   function will read the first numEntInInput entries of
///   lclColInds, sortPerm, and vals.
/// \param forceAtomic [in] Whether to use atomic updates when
///   modifying the entries of the matrix A.  For forwards
///   compatibility, this should always be a compile-time constant.
///   It defaults to whether the matrix's Kokkos execution space is
///   NOT Kokkos::Serial.
/// \param checkInputIndices [in] Whether to check whether the input
///   indices are valid column indices before just using them.  This
///   MUST be a compile-time constant.  Default is true, that is,
///   always check.
///
/// \return If checkInputIndices is true, return the number of input
///   indices that are valid column indices in that row of the matrix.
///   If checkInputIndices is false, just return numEntInInput.
template<class SparseMatrixType,
         class ValsViewType>
KOKKOS_FUNCTION
typename SparseMatrixType::ordinal_type
crsMatrixReplaceValues_sortedSortedLinear (const SparseMatrixType& A,
                                           const typename SparseMatrixType::ordinal_type lclRow,
                                           const typename SparseMatrixType::ordinal_type lclColInds[],
                                           const typename SparseMatrixType::ordinal_type sortPerm[],
                                           const ValsViewType& vals,
                                           const typename SparseMatrixType::ordinal_type numEntInInput,
                                           const bool forceAtomic =
#ifdef KOKKOS_ENABLE_SERIAL
                                           ! std::is_same<typename SparseMatrixType::device_type::execution_space, Kokkos::Serial>::type,
#else // NOT KOKKOS_ENABLE_SERIAL
                                           false,
#endif // KOKKOS_ENABLE_SERIAL
                                           const bool checkInputIndices = true)
{
  typedef typename std::remove_const<typename SparseMatrixType::value_type>::type
    matrix_scalar_type;
  static_assert (std::is_same<matrix_scalar_type,
                 typename SparseMatrixType::value_type>::value,
                 "The matrix's entries must have a nonconst type.");
  static_assert (std::is_assignable<matrix_scalar_type,
                 typename std::decay< decltype (A.values[0] + vals[0]) >::type>::value,
                 "The result of adding a matrix entry and an entry of vals "
                 "MUST be assignable to a matrix entry.");
  typedef typename SparseMatrixType::ordinal_type LO;
  static_assert (std::is_integral<LO>::value, "SparseMatrixType::ordinal_type "
                 "must be a built-in integer type.");

  // If lclRow is NOT a valid row index, this will return a view of
  // zero entries.  If checkInputIndices is true, thus, then none of
  // the input indices will be valid in that case.
  auto row_view = A.row (lclRow);
  const LO numEntInRow = static_cast<LO> (row_view.length);
  // Number of valid local column indices found, that is, the number
  // of input indices that are valid column indices found in row
  // lclRow of the matrix.  If not checking, we just return the number
  // of input indices.
  LO numValid = checkInputIndices ? static_cast<LO> (0) : numEntInRow;

  // Since both the matrix row and the input (after permutation) are
  // sorted, we only need to pass once over the matrix row.  'offset'
  // tells us the current search position in the matrix row.
  LO offset = 0;
  for (LO j = 0; j < numEntInInput; ++j) {
    const LO perm_index = sortPerm[j];
    const LO lclColInd = lclColInds[perm_index];
    // Search linearly in the matrix row for the current index.
    // If we ever want binary search, this would be the place.
    while (row_view.colidx(offset) != lclColInd) {
      ++offset;
    }

    // If checkInputIndices were a compile-time constant, then the
    // compiler might not need to insert a branch here.  This should
    // help vectorization, if vectorization is possible at all.
    if (checkInputIndices) {
      if (offset != numEntInRow) {
        // If forceAtomic were a compile-time constant, then the
        // compiler might not need to insert a branch here.  This
        // could help vectorization, if vectorization is possible.
        if (forceAtomic) {
          Kokkos::atomic_assign (&(row_view.value(offset)), vals[perm_index]);
        }
        else {
          row_view.value(offset) += vals[perm_index];
        }
        ++numValid;
      }
    }
    else { // don't check input indices; assume they are in the row
      // See above note on forceAtomic.
      if (forceAtomic) {
        Kokkos::atomic_add (&(row_view.value(offset)), vals[perm_index]);
      }
      else {
        row_view.value(offset) += vals[perm_index];
      }
    }
  }

  return numValid;
}

/// \brief <tt>A(lids[j], lids[j]) += lhs(j,j)</tt> and
///   <tt>x(lids[j]) += rhs(j)</tt>,
///   for all j in </tt>0 .. eltDim-1</tt>.
///
/// Assume the following:
/// <ul>
/// <li>In each row of the sparse matrix A, the column indices are
///   sorted.</li>
/// <li>The row and column indices of A have the same local indexing
///   scheme, that is, a valid row index is a valid column index and
///   vice versa.</li>
/// </ul>
/// Sum the dense "element" matrix (2-D Kokkos::View) \c lhs into the
/// entries of the sparse matrix A corresponding to the input row and
/// column indices \c lids.  Also, sum the dense "element" vector (1-D
/// Kokkos::View) \c rhs into the entries of the dense vector x
/// corresponding to the input row indices \c lids.
///
/// \tparam SparseMatrixType Specialization of KokkosSparse::CrsMatrix.
/// \tparam RhsViewType Specialization of a 1-D Kokkos::View.
/// \tparam LhsViewType Specialization of a 2-D Kokkos::View.
///
/// \param A [in/out] Sparse matrix (KokkosSparse::CrsMatrix) to modify.
/// \param x [in/out] Dense vector (1-D Kokkos::View) to modify.
/// \param lids [in/out] Local row and column indices of A to modify.
///   This function may sort this array, and output the permutation
///   that makes it sorted to \c sortPerm.  \c lids must have the same
///   number of entries as <tt>rhs.extent(0)</tt>,
///   <tt>lhs.extent(0)</tt>, and <tt>lhs.extent(1)</tt>.
/// \param sortPerm [out] Permutation that makes \c lids (on input)
///   sorted.  It must have the same number of writeable entries as
///   \c lids (see above).
/// \param rhs [in] Dense "element" vector of input values to sum into
///   the dense vector x; a 1-D Kokkos::View.  It must have the same
///   number of entries as each dimension of \c lhs.
/// \param lhs [in] Dense, square "element" matrix of input values to
///   sum into the sparse matrix A; a 2-D Kokkos::View.  Each of its
///   dimensions must be the same as the number of entries in \c rhs.
/// \param forceAtomic [in] Whether to use atomic updates when
///   modifying the entries of the matrix A and vector x.  For
///   forwards compatibility, this should always be a compile-time
///   constant.  It defaults to whether the matrix's Kokkos execution
///   space is NOT Kokkos::Serial.
/// \param checkInputIndices [in] Whether to check whether the input
///   indices are valid column indices before just using them.  This
///   MUST be a compile-time constant.  Default is true, that is,
///   always check.
///
/// \return If checkInputIndices is true, return the number of input
///   indices that are valid column indices in that row of the matrix.
///   If checkInputIndices is false, just return numEntInInput.
template<class SparseMatrixType,
         class VectorViewType,
         class RhsViewType,
         class LhsViewType>
KOKKOS_FUNCTION
typename SparseMatrixType::ordinal_type
crsMatrixAssembleElement_sortedLinear (const SparseMatrixType& A,
                                       const VectorViewType& x,
                                       typename SparseMatrixType::ordinal_type lids[],
                                       typename SparseMatrixType::ordinal_type sortPerm[],
                                       const RhsViewType& rhs,
                                       const LhsViewType& lhs,
                                       const bool forceAtomic =
#ifdef KOKKOS_ENABLE_SERIAL
                                       ! std::is_same<typename SparseMatrixType::device_type::execution_space, Kokkos::Serial>::type,
#else // NOT KOKKOS_ENABLE_SERIAL
                                       false,
#endif // KOKKOS_ENABLE_SERIAL
                                       const bool checkInputIndices = true)
{
  typedef typename std::remove_const<typename SparseMatrixType::value_type>::type
    matrix_scalar_type;
  typedef typename std::remove_const<typename VectorViewType::value_type>::type
    vector_scalar_type;
  static_assert (std::is_same<matrix_scalar_type,
                 typename SparseMatrixType::value_type>::value,
                 "The sparse output matrix A's entries must have a nonconst type.");
  static_assert (std::is_same<vector_scalar_type,
                 typename VectorViewType::value_type>::value,
                 "The dense output vector x's entries must have a nonconst type.");
  // static_assert (std::is_assignable<matrix_scalar_type,
  //                typename std::decay< decltype (A.values[0] + lhs(0,0)) >::type>::value,
  //                "The result of adding a sparse matrix entry and an entry of "
  //                "lhs (the dense element matrix) "
  //                "MUST be assignable to a matrix entry.");
  // static_assert (std::is_assignable<vector_scalar_type,
  //                typename std::decay< decltype (x[0] + rhs[0]) >::type>::value,
  //                "The result of adding a vector entry and an entry of "
  //                "rhs (the dense element vector) "
  //                "MUST be assignable to a vector entry.");
  typedef typename SparseMatrixType::ordinal_type LO;
  static_assert (std::is_integral<LO>::value, "SparseMatrixType::ordinal_type "
                 "must be a built-in integer type.");

  const LO eltDim = rhs.extent (0);

  // Generate sort permutation
  for (LO i = 0; i < eltDim; ++i) {
    sortPerm[i] = i;
  }
  shellSortKeysAndValues (lids, sortPerm, eltDim);

  LO totalNumValid = 0;
  for (LO r = 0; r < eltDim; ++r) {
    const LO lid = lids[r];
    //auto lhs_r = Kokkos::subview (lhs, sortPerm[r], Kokkos::ALL ());
    auto lhs_r = Kokkos::subview (lhs, r, Kokkos::ALL ());

    // This assumes that lid is always a valid row in the sparse
    // matrix, and that the local indices in each row of the matrix
    // are always sorted.
    const LO curNumValid =
      crsMatrixSumIntoValues_sortedSortedLinear (A, lid, lids, sortPerm, lhs_r,
                                                 eltDim, forceAtomic,
                                                 checkInputIndices);
    if (forceAtomic) {
      Kokkos::atomic_add (&x(lid), rhs(sortPerm[r]));
    }
    else {
      x(lid) += rhs(sortPerm[r]);
    }
    totalNumValid += curNumValid;
  }
  return totalNumValid;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSMATRIXASSEMBLEELEMENT_HPP
