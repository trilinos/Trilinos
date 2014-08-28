/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_SEQUENTIAL_SPARSEKERNELS_HPP
#define KOKKOS_SEQUENTIAL_SPARSEKERNELS_HPP

/// \file Kokkos_Sequential_SparseKernels.hpp
/// \brief Sequential implementations of (local) sparse kernels.
///
/// This file exists mainly as a temporary porting aid.  Until we can
/// set up thread-parallel versions of these kernels, we have
/// sequential versions here.  They are not terribly well optimized.
/// The point is to have reasonable sequential defaults, not
/// super-fast kernels.  If you want super-fast kernels, reimplement
/// them in Kokkos and try to make them parallel.
///
/// The exception to this might be Gauss-Seidel.  Gauss-Seidel is
/// harder to parallelize without changing the semantics, unless we
/// implement it with triangular solves and do the latter in parallel.
/// There is a set-up cost for thread-parallel sparse triangular
/// solves, and speed-up isn't even close to perfect, so it might pay
/// for smaller thread counts to have an optimized sequential kernel.
/// We have <i>not</i> done this here.

#include <KokkosLinAlg_config.h>
#include <Kokkos_ArithTraits.hpp>

namespace Kokkos {
namespace Sequential {

/// \brief Implementation of local Gauss-Seidel.
///
/// "Local" means local to the MPI process.
///
/// \tparam LocalOrdinal The type of each column index.
/// \tparam OffsetType The type of the entries of the row offsets array.
/// \tparam MatrixScalar The type of the entries (values of the matrix.
/// \tparam DomainScalar The type of the entries of the input multivector.
/// \tparam RangeScalar The type of the entries of the output multivector.
///
/// \param numRows [in] Number of rows in the (local) matrix.
/// \param numCols [in] Number of columns in the input multivector B.
///   <i>NOT</i> the number of columns in the matrix!
/// \param ptr [in] The matrix's row offsets.
/// \param ind [in] The matrix's column indices.
/// \param val [in] The matrix's values.
/// \param B [in] The input multivector.
/// \param b_stride [in] Column stride of the input multivector B.
/// \param X [in] The output multivector.
/// \param x_stride [in] Column stride of the output multivector X.
/// \param D [in] Array of the "diagonal entries" of the matrix.
///   These may differ from the matrix's actual diagonal entries, as
///   in L1 Gauss-Seidel for example.
/// \param omega [in] Damping parameter.
/// \param direction [in] Sweep direction.
///
template<class LocalOrdinal,
         class OffsetType,
         class MatrixScalar,
         class DomainScalar,
         class RangeScalar>
void
gaussSeidel (const LocalOrdinal numRows,
             const LocalOrdinal numCols,
             const OffsetType* const ptr,
             const LocalOrdinal* const ind,
             const MatrixScalar* const val,
             const DomainScalar* const B,
             const OffsetType b_stride,
             RangeScalar* const X,
             const OffsetType x_stride,
             const MatrixScalar* const D,
             const MatrixScalar omega,
             const Kokkos::ESweepDirection direction)
{
  using Kokkos::Details::ArithTraits;
  typedef LocalOrdinal LO;

  if (numRows == 0 || numCols == 0) {
    return; // Nothing to do.
  }
  else if (numRows > 0 && ptr[numRows] == 0) {
    // All the off-diagonal entries of A are zero, and all the
    // diagonal entries are (implicitly) 1.  Therefore compute: X :=
    // (1 - omega) X + omega B.  There's no need to care about the
    // direction, since there are no cross-row data dependencies in
    // this case.
    const MatrixScalar oneMinusOmega =
      ArithTraits<MatrixScalar>::one () - omega;
    for (OffsetType j = 0; j < numCols; ++j) {
      RangeScalar* const x_j = X + j*x_stride;
      const DomainScalar* const b_j = B + j*b_stride;
      for (OffsetType i = 0; i < numRows; ++i) {
        x_j[i] = oneMinusOmega * x_j[i] + omega * b_j[i];
      }
    }
    return;
  }

  if (numCols == 1) {
    if (direction == Kokkos::Forward) {
      for (LO i = 0; i < numRows; ++i) {
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    } else if (direction == Kokkos::Backward) {
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO i = numRows - 1; i != 0; --i) {
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
      { // last loop iteration
        const LO i = 0;
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    }
  }
  else { // numCols > 1
    // mfh 20 Dec 2012: If Gauss-Seidel for multivectors with
    // multiple columns becomes important, we can add unrolled
    // implementations.  The implementation below is not unrolled.
    // It may also be reasonable to parallelize over right-hand
    // sides, if there are enough of them, especially if the matrix
    // fits in cache.
    Teuchos::Array<RangeScalar> temp (numCols);
    RangeScalar* const x_temp = temp.getRawPtr ();

    if (direction == Kokkos::Forward) {
      for (LO i = 0; i < numRows; ++i) {
        for (OffsetType c = 0; c < numCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
    } else if (direction == Kokkos::Backward) { // backward mode
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO i = numRows - 1; i != 0; --i) {
        for (OffsetType c = 0; c < numCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
      { // last loop iteration
        const LO i = 0;
        for (OffsetType c = 0; c < numCols; ++c) {
          x_temp[c] = ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
    }
  }
}


/// \brief Implementation of reordered local Gauss-Seidel.
///
/// "Local" means local to the MPI process.
///
/// \tparam LocalOrdinal The type of each column index.
/// \tparam OffsetType The type of the entries of the row offsets array.
/// \tparam MatrixScalar The type of the entries (values of the matrix.
/// \tparam DomainScalar The type of the entries of the input multivector.
/// \tparam RangeScalar The type of the entries of the output multivector.
///
/// \param numRows [in] Number of rows in the (local) matrix.
/// \param numCols [in] Number of columns in the input multivector B.
///   <i>NOT</i> the number of columns in the matrix!
/// \param ptr [in] The matrix's row offsets.
/// \param ind [in] The matrix's column indices.
/// \param val [in] The matrix's values.
/// \param B [in] The input multivector.
/// \param b_stride [in] Column stride of the input multivector B.
/// \param X [in] The output multivector.
/// \param x_stride [in] Column stride of the output multivector X.
/// \param D [in] Array of the "diagonal entries" of the matrix.
///   These may differ from the matrix's actual diagonal entries, as
///   in L1 Gauss-Seidel for example.
/// \param rowInd [in] Array of row indices to process.  It has
///   numRowInds entries.  This array determines the order in which
///   the rows are accessed.  It is legal for this to contain fewer
///   entries than the number of rows in the matrix.
/// \param numRowInds [in] Number of entries in rowInd; the number of
///   rows to process.  This may be less than or equal to numRows (the
///   number of rows in the matrix).
/// \param omega [in] Damping parameter.
/// \param direction [in] Sweep direction.
///
template<class LocalOrdinal,
         class OffsetType,
         class MatrixScalar,
         class DomainScalar,
         class RangeScalar>
void
reorderedGaussSeidel (const LocalOrdinal numRows,
                      const LocalOrdinal numCols,
                      const OffsetType* const ptr,
                      const LocalOrdinal* const ind,
                      const MatrixScalar* const val,
                      const DomainScalar* const B,
                      const OffsetType b_stride,
                      RangeScalar* const X,
                      const OffsetType x_stride,
                      const MatrixScalar* const D,
                      const LocalOrdinal* const rowInd,
                      const LocalOrdinal numRowInds, // length of rowInd
                      const MatrixScalar omega,
                      const Kokkos::ESweepDirection direction)
{
  using Kokkos::Details::ArithTraits;
  typedef LocalOrdinal LO;

  if (numRows == 0 || numCols == 0) {
    return; // Nothing to do.
  }
  else if (numRows > 0 && ptr[numRows] == 0) {
    // All the off-diagonal entries of A are zero, and all the
    // diagonal entries are (implicitly) 1.  Therefore compute: X :=
    // (1 - omega) X + omega B.  There's no need to care about the
    // direction or row ordering, since there are no cross-row data
    // dependencies in this case.
    const MatrixScalar oneMinusOmega =
      ArithTraits<MatrixScalar>::one () - omega;
    for (OffsetType j = 0; j < numCols; ++j) {
      RangeScalar* const x_j = X + j*x_stride;
      const DomainScalar* const b_j = B + j*b_stride;
      for (OffsetType i = 0; i < numRows; ++i) {
        x_j[i] = oneMinusOmega * x_j[i] + omega * b_j[i];
      }
    }
    return;
  }

  if (numCols == 1) {
    if (direction == Kokkos::Forward) {
      for (LO ii = 0; ii < numRowInds; ++ii) {
        LO i = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    } else if (direction == Kokkos::Backward) {
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO ii = numRowInds - 1; ii != 0; --ii) {
        LO i = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
      { // last loop iteration
        const LO ii = 0;
        LO i = rowInd[ii];
        RangeScalar x_temp = ArithTraits<RangeScalar>::zero ();
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          x_temp += A_ij * X[j];
        }
        X[i] += omega * D[i] * (B[i] - x_temp);
      }
    }
  }
  else { // numCols > 1
    // mfh 20 Dec 2012: If Gauss-Seidel for multivectors with
    // multiple columns becomes important, we can add unrolled
    // implementations.  The implementation below is not unrolled.
    // It may also be reasonable to parallelize over right-hand
    // sides, if there are enough of them, especially if the matrix
    // fits in cache.
    Teuchos::Array<RangeScalar> temp (numCols);
    RangeScalar* const x_temp = temp.getRawPtr ();

    if (direction == Kokkos::Forward) {
      for (LO ii = 0; ii < numRowInds; ++ii) {
        LO i = rowInd[ii];
        for (size_t c = 0; c < numCols; ++c) {
          x_temp[c] = Kokkos::Details::ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
    } else if (direction == Kokkos::Backward) { // backward mode
      // Split the loop so that it is correct even if LO is unsigned.
      // It's a bad idea for LO to be unsigned, but we want this to
      // work nevertheless.
      for (LO ii = numRowInds - 1; ii != 0; --ii) {
        LO i = rowInd[ii];
        for (size_t c = 0; c < numCols; ++c) {
          x_temp[c] = Kokkos::Details::ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
      { // last loop iteration
        const LO ii = 0;
        LO i = rowInd[ii];
        for (size_t c = 0; c < numCols; ++c) {
          x_temp[c] = Kokkos::Details::ArithTraits<RangeScalar>::zero ();
        }
        for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
          const LO j = ind[k];
          const MatrixScalar A_ij = val[k];
          for (OffsetType c = 0; c < numCols; ++c) {
            x_temp[c] += A_ij * X[j + x_stride*c];
          }
        }
        for (OffsetType c = 0; c < numCols; ++c) {
          X[i + x_stride*c] += omega * D[i] * (B[i + b_stride*c] - x_temp[c]);
        }
      }
    }
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsrUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type r = 0; r < numRows; ++r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // for each row r
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsr (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type r = 0; r < numRows; ++r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }

    matrix_scalar_type A_rr = STS::zero ();
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);

    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  The
      // way to fix this would be for Tpetra::CrsMatrix to remember
      // the local column index of the diagonal entry (if there is
      // one) in each row, and pass that along to this function.
      if (r == c) {
        A_rr += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // for each row r
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsrUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  // If local_ordinal_type is unsigned and numRows is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numRows == 0) {
    return;
  }

  // Don't use r >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do r == 0 (last
  // iteration) below.
  for (local_ordinal_type r = numRows - 1; r != 0; --r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // for each row r

  // Last iteration: r = 0.
  {
    const local_ordinal_type r = 0;
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
  } // last iteration: r = 0
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsr (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  //const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  // If local_ordinal_type is unsigned and numRows is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numRows == 0) {
    return;
  }

  // Don't use r >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do r == 0 (last
  // iteration) below.
  for (local_ordinal_type r = numRows - 1; r != 0; --r) {
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    // We assume the diagonal entry is first in the row.
    const matrix_scalar_type A_rr = val(beg);
    for (offset_type k = beg + static_cast<offset_type> (1); k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // for each row r

  // Last iteration: r = 0.
  {
    const local_ordinal_type r = 0;
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) = Y(r, j);
    }
    const offset_type beg = ptr(r);
    const offset_type end = ptr(r+1);
    // We assume the diagonal entry is first in the row.
    const matrix_scalar_type A_rr = val(beg);
    for (size_t k = beg + 1; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type c = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current row r
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(r, j) /= A_rr;
    }
  } // last iteration: r = 0
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type r = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const matrix_scalar_type A_rc = val(k);
      const local_ordinal_type r = ind(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCsc (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscUnitDiag (RangeMultiVectorType X,
                          const CrsMatrixType& A,
                          DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscUnitDiagConj (RangeMultiVectorType X,
                              const CrsMatrixType& A,
                              DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
upperTriSolveCscConj (RangeMultiVectorType X,
                      const CrsMatrixType& A,
                      DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  // If local_ordinal_type is unsigned and numCols is 0, the loop
  // below will have entirely the wrong number of iterations.
  if (numCols == 0) {
    return;
  }

  // Don't use c >= 0 as the test, because that fails if
  // local_ordinal_type is unsigned.  We do c == 0 (last
  // iteration) below.
  for (local_ordinal_type c = numCols - 1; c != 0; --c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c

  // Last iteration: c = 0.
  {
    const local_ordinal_type c = 0;
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  }
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCsc (RangeMultiVectorType X,
                  const CrsMatrixType& A,
                  DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = val(k);
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscUnitDiagConj (RangeMultiVectorType X,
                              const CrsMatrixType& A,
                              DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      for (local_ordinal_type j = 0; j < numVecs; ++j) {
        X(r, j) -= A_rc * X(c, j);
      }
    } // for each entry A_rc in the current column c
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
lowerTriSolveCscConj (RangeMultiVectorType X,
                      const CrsMatrixType& A,
                      DomainMultiVectorType Y)
{
  typedef typename CrsMatrixType::row_map_type::non_const_value_type offset_type;
  typedef typename CrsMatrixType::index_type::non_const_value_type local_ordinal_type;
  typedef typename CrsMatrixType::values_type::non_const_value_type matrix_scalar_type;
  typedef Kokkos::Details::ArithTraits<matrix_scalar_type> STS;

  const local_ordinal_type numRows = A.numRows ();
  const local_ordinal_type numCols = A.numCols ();
  const local_ordinal_type numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;
  typename CrsMatrixType::index_type ind = A.graph.entries;
  typename CrsMatrixType::values_type val = A.values;

  for (local_ordinal_type j = 0; j < numVecs; ++j) {
    for (local_ordinal_type i = 0; i < numRows; ++i) {
      X(i, j) = Y(i, j);
    }
  }

  for (local_ordinal_type c = 0; c < numCols; ++c) {
    matrix_scalar_type A_cc = STS::zero ();
    const offset_type beg = ptr(c);
    const offset_type end = ptr(c+1);
    for (offset_type k = beg; k < end; ++k) {
      const local_ordinal_type r = ind(k);
      const matrix_scalar_type A_rc = STS::conj (val(k));
      // FIXME (mfh 28 Aug 2014) This assumes that the diagonal entry
      // has equal local row and column indices.  That may not
      // necessarily hold, depending on the row and column Maps.  See
      // note above.
      if (r == c) {
        A_cc += A_rc;
      } else {
        for (local_ordinal_type j = 0; j < numVecs; ++j) {
          X(r, j) -= A_rc * X(c, j);
        }
      }
    } // for each entry A_rc in the current column c
    for (local_ordinal_type j = 0; j < numVecs; ++j) {
      X(c, j) /= A_cc;
    }
  } // for each column c
}


template<class CrsMatrixType,
         class DomainMultiVectorType,
         class RangeMultiVectorType>
void
triSolveKokkos (RangeMultiVectorType X,
                const CrsMatrixType& A,
                DomainMultiVectorType Y,
                const Teuchos::EUplo triUplo,
                const Teuchos::EDiag unitDiag,
                const Teuchos::ETransp trans)
{
  typedef typename CrsMatrixType::index_type::non_const_value_type LO;
  const char prefix[] = "Kokkos::Sequential::triSolveKokkos: ";
  const LO numRows = A.numRows ();
  const LO numCols = A.numCols ();
  const LO numVecs = X.dimension_1 ();
  typename CrsMatrixType::row_map_type ptr = A.graph.row_map;

  TEUCHOS_TEST_FOR_EXCEPTION(
    triUplo != Teuchos::LOWER_TRI && triUplo != Teuchos::UPPER_TRI &&
    triUplo != Teuchos::UNDEF_TRI,
    std::invalid_argument, prefix << "triUplo has an invalid value " << triUplo
    << ".  Valid values are Teuchos::LOWER_TRI=" << Teuchos::LOWER_TRI <<
    ", Teuchos::UPPER_TRI=" << Teuchos::UPPER_TRI << ", and Teuchos::UNDEF_TRI="
    << Teuchos::UNDEF_TRI << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    triUplo == Teuchos::UNDEF_TRI, std::invalid_argument, prefix <<
    "The matrix is neither lower nor upper triangular (triUplo="
    "Teuchos::UNDEF_TRI), so you may not call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    unitDiag != Teuchos::UNIT_DIAG && unitDiag != Teuchos::NON_UNIT_DIAG,
    std::invalid_argument, prefix << "unitDiag has an invalid value "
    << unitDiag << ".  Valid values are Teuchos::UNIT_DIAG="
    << Teuchos::UNIT_DIAG << " and Teuchos::NON_UNIT_DIAG="
    << Teuchos::NON_UNIT_DIAG << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    unitDiag != Teuchos::UNIT_DIAG && numRows > 0 && ptr(numRows) == 0,
    std::invalid_argument, prefix << "Triangular solve with an empty matrix "
    "is only valid if the matrix has an implicit unit diagonal.  This matrix "
    "does not.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    trans != Teuchos::NO_TRANS && trans != Teuchos::TRANS &&
    trans != Teuchos::CONJ_TRANS,
    std::invalid_argument, prefix << "trans has an invalid value " << trans
    << ".  Valid values are Teuchos::NO_TRANS=" << Teuchos::NO_TRANS << ", "
    << "Teuchos::TRANS=" << Teuchos::TRANS << ", and Teuchos::CONJ_TRANS="
    << Teuchos::CONJ_TRANS << ".");

  TEUCHOS_TEST_FOR_EXCEPTION(
    numRows != X.dimension_0 (), std::invalid_argument, prefix << "numRows = "
    << numRows << " != X.dimension_0() = " << X.dimension_0 () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numCols != Y.dimension_0 (), std::invalid_argument, prefix << "numCols = "
    << numCols << " != Y.dimension_0() = " << Y.dimension_0 () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.dimension_1 (), std::invalid_argument, prefix <<
    "X.dimension_1 () = " << numVecs << " != Y.dimension_1 () = "
    << Y.dimension_1 () << ".");

  if (trans == Teuchos::NO_TRANS) {          // no transpose
    if (triUplo == Teuchos::LOWER_TRI) { // lower triangular
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        lowerTriSolveCsrUnitDiag (X, A, Y);
      } else {                          // non unit diagonal
        lowerTriSolveCsr (X, A, Y);
      }
    } else {                             // upper triangular
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        upperTriSolveCsrUnitDiag (X, A, Y);
      } else {                          // non unit diagonal
        upperTriSolveCsr (X, A, Y);
      }
    }
  }
  else if (trans == Teuchos::TRANS) {           // transpose
    if (triUplo == Teuchos::LOWER_TRI) { // lower triangular
      // Transposed lower tri CSR => upper tri CSC.
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        upperTriSolveCscUnitDiag (X, A, Y);
      } else {                          // non unit diagonal
        upperTriSolveCsc (X, A, Y);
      }
    }
    else {                               // upper triangular
      // Transposed upper tri CSR => lower tri CSC.
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        lowerTriSolveCscUnitDiag (X, A, Y);
      } else {                          // non unit diagonal
        lowerTriSolveCsc (X, A, Y);
      }
    }
  }
  else if (trans == Teuchos::CONJ_TRANS) { // conj transpose
    if (triUplo == Teuchos::LOWER_TRI) { // lower triangular
      // Transposed lower tri CSR => upper tri CSC.
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        upperTriSolveCscUnitDiagConj (X, A, Y);
      } else {                          // non unit diagonal
        upperTriSolveCscConj (X, A, Y);
      }
    }
    else {                               // upper triangular
      // Transposed upper tri CSR => lower tri CSC.
      if (unitDiag == Teuchos::UNIT_DIAG) { // unit diagonal
        lowerTriSolveCscUnitDiagConj (X, A, Y);
      } else {                          // non unit diagonal
        lowerTriSolveCscConj (X, A, Y);
      }
    }
  }
}


} // namespace Sequential
} // namespace Kokkos

#endif // KOKKOS_SEQUENTIAL_SPARSEKERNELS_HPP
